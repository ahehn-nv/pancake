// Authors: Ivan Sovic

#include <pacbio/alignment/AlignmentTools.h>
#include <pacbio/pancake/AlignerBase.h>
#include <pacbio/pancake/MapperBatch.h>
#include <pacbio/pancake/OverlapWriterFactory.h>
#include <pacbio/util/TicToc.h>
#include <pbcopper/logging/Logging.h>
#include <pbcopper/parallel/FireAndForget.h>
#include <pbcopper/parallel/WorkQueue.h>
#include <algorithm>
#include <array>
#include <iostream>
#include <sstream>
#include <tuple>

namespace PacBio {
namespace Pancake {

MapperBatch::MapperBatch(const MapperCLRSettings& settings, int32_t numThreads)
    : settings_{settings}, numThreads_(numThreads)
{
}

MapperBatch::~MapperBatch() = default;

std::vector<std::vector<MapperBaseResult>> MapperBatch::DummyMapAndAlign(
    const std::vector<MapperBatchChunk>& batchData)
{
    return DummyMapAndAlignImpl_(batchData, settings_, numThreads_);
}

std::vector<std::vector<MapperBaseResult>> MapperBatch::MapAndAlignCPU(
    const std::vector<MapperBatchChunk>& batchData)
{
    return MapAndAlignCPUImpl_(batchData, settings_, numThreads_);
}

std::vector<std::vector<MapperBaseResult>> MapperBatch::DummyMapAndAlignImpl_(
    const std::vector<MapperBatchChunk>& batchChunks, MapperCLRSettings settings,
    int32_t /*numThreads*/)
{
    std::vector<std::vector<MapperBaseResult>> results;
    results.reserve(batchChunks.size());

    // Deactivate alignment for the mapper.
    MapperCLRSettings settingsCopy = settings;
    settingsCopy.align = false;
    auto mapper = std::make_unique<MapperCLR>(settingsCopy);

    // Map.
    for (size_t i = 0; i < batchChunks.size(); ++i) {
        const auto& chunk = batchChunks[i];
        std::vector<MapperBaseResult> result =
            mapper->MapAndAlign(chunk.targetSeqs, chunk.querySeqs);
        results.emplace_back(std::move(result));
    }

    if (settings.align) {
        for (size_t i = 0; i < results.size(); ++i) {
            auto& result = results[i];
            const auto& chunk = batchChunks[i];
            for (size_t qId = 0; qId < result.size(); ++qId) {
                result[qId] = mapper->Align(chunk.targetSeqs, chunk.querySeqs[qId], result[qId]);
            }
        }
    }

    return results;
}

std::vector<std::vector<MapperBaseResult>> MapperBatch::MapAndAlignCPUImpl_(
    const std::vector<MapperBatchChunk>& batchChunks, MapperCLRSettings settings,
    int32_t numThreads)
{
    std::vector<std::vector<MapperBaseResult>> results(batchChunks.size());

    // Deactivate alignment for the mapper.
    MapperCLRSettings settingsCopy = settings;
    settingsCopy.align = false;

    // Determine how many records should land in each thread, spread roughly evenly.
    const int32_t numRecords = batchChunks.size();
    const int32_t actualThreadCount = std::min(static_cast<int32_t>(numThreads), numRecords);
    const int32_t minimumRecordsPerThreads = (numRecords / actualThreadCount);
    const int32_t remainingRecords = (numRecords % actualThreadCount);
    std::vector<int32_t> recordsPerThread(actualThreadCount, minimumRecordsPerThreads);
    for (int32_t i = 0; i < remainingRecords; ++i) {
        ++recordsPerThread[i];
    }

    // Run the mapping in parallel.
    PacBio::Parallel::FireAndForget faf(numThreads);
    int32_t submittedCount = 0;
    for (int32_t i = 0; i < actualThreadCount; ++i) {
        faf.ProduceWith(WorkerMapper_, std::cref(batchChunks), std::cref(settingsCopy),
                        submittedCount, submittedCount + recordsPerThread[i], std::ref(results));
        submittedCount += recordsPerThread[i];
    }
    faf.Finalize();

    if (settings.align) {
        AlignerBatchCPU alignerInternal(settings.alignerTypeGlobal, settings.alnParamsGlobal,
                                        settings.alignerTypeExt, settings.alnParamsExt);
        PrepareSequencesForAlignment_(alignerInternal, batchChunks, results,
                                      BatchAlignerRegionType::GLOBAL);

        AlignerBatchCPU alignerFlanks(settings.alignerTypeGlobal, settings.alnParamsGlobal,
                                      settings.alignerTypeExt, settings.alnParamsExt);
        PrepareSequencesForAlignment_(alignerFlanks, batchChunks, results,
                                      BatchAlignerRegionType::SEMIGLOBAL);

        alignerInternal.AlignAll(numThreads);
        alignerFlanks.AlignAll(numThreads);

        StitchAlignments_(batchChunks, results, alignerInternal.GetAlnResults(),
                          alignerFlanks.GetAlnResults());

        UpdateSecondaryAndFilter_(results, settings.secondaryAllowedOverlapFractionQuery,
                                  settings.secondaryAllowedOverlapFractionTarget,
                                  settings.secondaryMinScoreFraction, settings.bestNSecondary);
    }

    return results;
}

void MapperBatch::WorkerMapper_(const std::vector<MapperBatchChunk>& batchChunks,
                                const MapperCLRSettings& settings, int32_t startId, int32_t endId,
                                std::vector<std::vector<MapperBaseResult>>& results)
{
    auto mapper = std::make_unique<MapperCLR>(settings);
    for (int32_t i = startId; i < endId; ++i) {
        const auto& chunk = batchChunks[i];
        results[i] = mapper->MapAndAlign(chunk.targetSeqs, chunk.querySeqs);
    }
}

void MapperBatch::PrepareSequencesForAlignment_(
    AlignerBatchCPU& aligner, const std::vector<MapperBatchChunk>& batchChunks,
    const std::vector<std::vector<MapperBaseResult>>& mappingResults,
    const BatchAlignerRegionType& regionsToAdd)
{
    // Results are a vector for every chunk (one chunk is one ZMW).
    for (size_t resultId = 0; resultId < mappingResults.size(); ++resultId) {
        const auto& result = mappingResults[resultId];
        const auto& chunk = batchChunks[resultId];
        // std::cerr << "[resultId = " << resultId << "]\n";

        // One chunk can have multiple queries (subreads).
        for (size_t qId = 0; qId < result.size(); ++qId) {
            // std::cerr << "    [qId = " << qId << "]\n";

            // Prepare the query data in fwd and rev.
            const char* qSeqFwd = chunk.querySeqs[qId].c_str();
            const int32_t qLen = chunk.querySeqs[qId].size();
            const std::string qSeqRevString = PacBio::Pancake::ReverseComplement(qSeqFwd, 0, qLen);
            const char* qSeqRev = qSeqRevString.c_str();

            // Each query can have multiple mappings.
            for (size_t mapId = 0; mapId < result[qId].mappings.size(); ++mapId) {
                const auto& mapping = result[qId].mappings[mapId];

                // std::cerr << "        [mapId = " << mapId << "] "
                //           << OverlapWriterBase::PrintOverlapAsM4(mapping->mapping, "", "", true,
                //                                                  true)
                //           << "\n";

                const char* tSeq = chunk.targetSeqs[mapping->mapping->Bid].c_str();
                const int32_t tLen = chunk.targetSeqs[mapping->mapping->Bid].size();

                // Each mapping is split into regions in between seed hits for alignment.
                std::string qSubSeq;
                std::string tSubSeq;
                for (size_t regId = 0; regId < mapping->regionsForAln.size(); ++regId) {
                    const auto& region = mapping->regionsForAln[regId];
                    // std::cerr << "            [regId = " << regId << "] " << region << "\n";

                    if ((region.type == RegionType::GLOBAL &&
                         regionsToAdd == BatchAlignerRegionType::SEMIGLOBAL) ||
                        (region.type != RegionType::GLOBAL &&
                         regionsToAdd == BatchAlignerRegionType::GLOBAL)) {
                        continue;
                    }

                    // Prepare the sequences for alignment.
                    const char* qSeqInStrand = region.queryRev ? qSeqRev : qSeqFwd;
                    const char* tSeqInStrand = tSeq;
                    int32_t qStart = region.qStart;
                    int32_t tStart = region.tStart;
                    const int32_t qSpan = region.qSpan;
                    const int32_t tSpan = region.tSpan;
                    if (region.type == RegionType::FRONT) {
                        qSubSeq = std::string(qSeqInStrand + qStart, qSpan);
                        tSubSeq = std::string(tSeqInStrand + tStart, tSpan);
                        std::reverse(qSubSeq.begin(), qSubSeq.end());
                        std::reverse(tSubSeq.begin(), tSubSeq.end());
                        qSeqInStrand = qSubSeq.c_str();
                        tSeqInStrand = tSubSeq.c_str();
                        qStart = 0;
                        tStart = 0;
                    }

                    // Add the region.
                    aligner.AddSequencePair(qSeqInStrand + qStart, qSpan, tSeqInStrand + tStart,
                                            tSpan, region.type == RegionType::GLOBAL);
                }
            }
        }
    }
}

void MapperBatch::StitchAlignments_(
    const std::vector<MapperBatchChunk>& batchChunks,
    const std::vector<std::vector<MapperBaseResult>>& mappingResults,
    const std::vector<AlignmentResult>& internalAlns, const std::vector<AlignmentResult>& flankAlns)
{
    size_t currInternal = 0;
    size_t currFlank = 0;

    // Results are a vector for every chunk (one chunk is one ZMW).
    for (size_t resultId = 0; resultId < mappingResults.size(); ++resultId) {
        auto& result = mappingResults[resultId];
        const auto& chunk = batchChunks[resultId];
        // std::cerr << "[stitch resultId = " << resultId << "]\n";

        // std::vector<PacBio::Pancake::MapperBaseResult>
        // One chunk can have multiple queries (subreads).
        for (size_t qId = 0; qId < result.size(); ++qId) {
            // std::cerr << "    [stitch qId = " << qId << "]\n";

            // Prepare the query data in fwd and rev.
            const char* qSeqFwd = chunk.querySeqs[qId].c_str();
            const int32_t qLen = chunk.querySeqs[qId].size();
            const std::string qSeqRevString = PacBio::Pancake::ReverseComplement(qSeqFwd, 0, qLen);
            const char* qSeqRev = qSeqRevString.c_str();

            // Each query can have multiple mappings.
            // std::vector<OverlapPtr> newMappings;
            for (size_t mapId = 0; mapId < result[qId].mappings.size(); ++mapId) {
                auto& mapping = result[qId].mappings[mapId];
                auto& ovl = mapping->mapping;
                ovl->Cigar.clear();
                int32_t newQueryStart = 0;
                int32_t newTargetStart = 0;
                int32_t newQueryEnd = 0;
                int32_t newTargetEnd = 0;

                // std::cerr << "        [stitch mapId = " << mapId << "] "
                //           << OverlapWriterBase::PrintOverlapAsM4(mapping->mapping, "", "", true,
                //                                                  true)
                //           << "\n";

                // Each mapping is split into regions in between seed hits for alignment.
                for (size_t regId = 0; regId < mapping->regionsForAln.size(); ++regId) {
                    const auto& region = mapping->regionsForAln[regId];
                    // std::cerr << "            [stitch regId = " << regId << "] " << region << "\n";

                    if (region.type == RegionType::FRONT) {
                        if (currFlank >= flankAlns.size()) {
                            std::ostringstream oss;
                            oss << "Invalid number of flank alignments. About to access currFlank "
                                   "= "
                                << currFlank << ", but flankAlns.size() = " << flankAlns.size()
                                << ". resultId = " << resultId << ", qId = " << qId
                                << ", mapId = " << mapId << ", regId = " << regId;
                            throw std::runtime_error(oss.str());
                        }
                        const auto& aln = flankAlns[currFlank];
                        // std::cerr << "                (front) [currFlank = " << currFlank
                        //           << ", currInternal = " << currInternal << "] aln = " << aln
                        //           << "\n";
                        size_t currCigarLen = ovl->Cigar.size();
                        ovl->Cigar.insert(ovl->Cigar.end(), aln.cigar.begin(), aln.cigar.end());
                        std::reverse(ovl->Cigar.begin() + currCigarLen, ovl->Cigar.end());
                        newQueryStart = region.qStart + region.qSpan - aln.lastQueryPos;
                        newTargetStart = region.tStart + region.tSpan - aln.lastTargetPos;
                        ++currFlank;
                    } else if (region.type == RegionType::BACK) {
                        if (currFlank >= flankAlns.size()) {
                            std::ostringstream oss;
                            oss << "Invalid number of flank alignments. About to access currFlank "
                                   "= "
                                << currFlank << ", but flankAlns.size() = " << flankAlns.size()
                                << ". resultId = " << resultId << ", qId = " << qId
                                << ", mapId = " << mapId << ", regId = " << regId;
                            throw std::runtime_error(oss.str());
                        }
                        const auto& aln = flankAlns[currFlank];
                        // std::cerr << "                (back) [currFlank = " << currFlank
                        //           << ", currInternal = " << currInternal << "] aln = " << aln
                        //           << "\n";
                        ovl->Cigar.insert(ovl->Cigar.end(), aln.cigar.begin(), aln.cigar.end());
                        newQueryEnd = region.qStart + aln.lastQueryPos;
                        newTargetEnd = region.tStart + aln.lastTargetPos;
                        ++currFlank;
                    } else {
                        if (currInternal >= internalAlns.size()) {
                            std::ostringstream oss;
                            oss << "Invalid number of internal alignments. About to access "
                                   "currInternal = "
                                << currInternal
                                << ", but internalAlns.size() = " << internalAlns.size()
                                << ". resultId = " << resultId << ", qId = " << qId
                                << ", mapId = " << mapId << ", regId = " << regId;
                            throw std::runtime_error(oss.str());
                        }
                        const auto& aln = internalAlns[currInternal];
                        // std::cerr << "                (internal) [currFlank = " << currFlank
                        //           << ", currInternal = " << currInternal << "] aln = " << aln
                        //           << "\n";
                        ovl->Cigar.insert(ovl->Cigar.end(), aln.cigar.begin(), aln.cigar.end());
                        ++currInternal;
                    }
                    // std::cerr << "CIGAR: " << ovl->Cigar.ToStdString() << "\n";
                }
                ovl->Astart = newQueryStart;
                ovl->Aend = newQueryEnd;
                ovl->Bstart = newTargetStart;
                ovl->Bend = newTargetEnd;

                // Reverse the CIGAR and the coordinates if needed.
                if (ovl->Brev) {
                    // CIGAR reversal.
                    std::reverse(ovl->Cigar.begin(), ovl->Cigar.end());

                    // Reverse the query coordinates.
                    std::swap(ovl->Astart, ovl->Aend);
                    ovl->Astart = ovl->Alen - ovl->Astart;
                    ovl->Aend = ovl->Alen - ovl->Aend;

                    // Get the forward-oriented target coordinates.
                    std::swap(ovl->Bstart, ovl->Bend);
                    ovl->Bstart = ovl->Blen - ovl->Bstart;
                    ovl->Bend = ovl->Blen - ovl->Bend;
                }

                // Set the alignment identity and edit distance.
                Alignment::DiffCounts diffs = CigarDiffCounts(ovl->Cigar);
                diffs.Identity(false, false, ovl->Identity, ovl->EditDistance);
                ovl->Score = -diffs.numEq;
                // std::cerr << "        [stitch after mapId = " << mapId << "] "
                //           << OverlapWriterBase::PrintOverlapAsM4(mapping->mapping, "", "", true,
                //                                                  true)
                //           << "\n";
            }
        }
    }
}

void MapperBatch::UpdateSecondaryAndFilter_(
    std::vector<std::vector<MapperBaseResult>>& mappingResults,
    double secondaryAllowedOverlapFractionQuery, double secondaryAllowedOverlapFractionTarget,
    double secondaryMinScoreFraction, int32_t bestNSecondary)
{
    // Results are a vector for every chunk (one chunk is one ZMW).
    for (size_t resultId = 0; resultId < mappingResults.size(); ++resultId) {
        auto& result = mappingResults[resultId];
        // One chunk can have multiple queries (subreads).
        for (size_t qId = 0; qId < result.size(); ++qId) {
            // Secondary/supplementary flagging.
            WrapFlagSecondaryAndSupplementary(
                result[qId].mappings, secondaryAllowedOverlapFractionQuery,
                secondaryAllowedOverlapFractionTarget, secondaryMinScoreFraction);
            CondenseMappings(result[qId].mappings, bestNSecondary);
        }
    }
}

}  // namespace Pancake
}  // namespace PacBio
