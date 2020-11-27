// Authors: Ivan Sovic

#include <pacbio/pancake/AlignerBase.h>
#include <pacbio/pancake/MapperBatch.h>
#include <pacbio/pancake/OverlapWriterFactory.h>
#include <pacbio/util/TicToc.h>
#include <pbcopper/logging/Logging.h>
#include <algorithm>
#include <array>
#include <iostream>
#include <sstream>
#include <tuple>

namespace PacBio {
namespace Pancake {

MapperBatch::MapperBatch(const MapperCLRSettings& settings, int32_t numThreads)
    : settings_{settings}, numThreads_(numThreads), mapper_(nullptr)
{
    // Deactivate alignment for the mapper.
    MapperCLRSettings settingsCopy = settings;
    settingsCopy.align = false;
    mapper_ = std::make_unique<MapperCLR>(settingsCopy);
}

MapperBatch::~MapperBatch() = default;

std::vector<std::vector<MapperBaseResult>> MapperBatch::DummyMapAndAlign(
    const std::vector<MapperBatchChunk>& batchData)
{
    return DummyMapAndAlignImpl_(mapper_, batchData, settings_, numThreads_);
}

std::vector<std::vector<MapperBaseResult>> MapperBatch::MapAndAlignCPU(
    const std::vector<MapperBatchChunk>& batchData)
{
    return MapAndAlignCPUImpl_(mapper_, batchData, settings_, numThreads_);
}

std::vector<std::vector<MapperBaseResult>> MapperBatch::DummyMapAndAlignImpl_(
    std::unique_ptr<MapperCLR>& mapper, const std::vector<MapperBatchChunk>& batchChunks,
    MapperCLRSettings settings, int32_t /*numThreads*/)
{
    std::vector<std::vector<MapperBaseResult>> results;
    results.reserve(batchChunks.size());

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
    std::unique_ptr<MapperCLR>& mapper, const std::vector<MapperBatchChunk>& batchChunks,
    MapperCLRSettings settings, int32_t numThreads)
{
    std::vector<std::vector<MapperBaseResult>> results;
    results.reserve(batchChunks.size());

    for (size_t i = 0; i < batchChunks.size(); ++i) {
        const auto& chunk = batchChunks[i];
        std::vector<MapperBaseResult> result =
            mapper->MapAndAlign(chunk.targetSeqs, chunk.querySeqs);
        results.emplace_back(std::move(result));
    }

    if (settings.align) {
        AlignerBasePtr alignerGlobal =
            AlignerFactory(settings.alignerTypeGlobal, settings.alnParamsGlobal);
        AlignerBasePtr alignerExt = AlignerFactory(settings.alignerTypeExt, settings.alnParamsExt);

        AlignerBatchCPU alignerInternal(alignerGlobal, alignerExt);
        PrepareSequencesForAlignment_(alignerInternal, batchChunks, results,
                                      BatchAlignerRegionType::GLOBAL);

        AlignerBatchCPU alignerFlanks(alignerGlobal, alignerExt);
        PrepareSequencesForAlignment_(alignerFlanks, batchChunks, results,
                                      BatchAlignerRegionType::SEMIGLOBAL);

        alignerInternal.AlignAll();
        alignerFlanks.AlignAll();

        results = StitchAlignments_(batchChunks, results, alignerInternal.GetAlnResults(),
                                    alignerFlanks.GetAlnResults());
    }

    return results;
}

void MapperBatch::PrepareSequencesForAlignment_(
    AlignerBatchCPU& aligner, const std::vector<MapperBatchChunk>& batchChunks,
    const std::vector<std::vector<MapperBaseResult>>& mappingResults,
    const BatchAlignerRegionType& regionsToAdd)
{
    // Results are a vector for every chunk (one chunk is one ZMW).
    for (size_t resultId = 0; resultId < mappingResults.size(); ++resultId) {
        auto& result = mappingResults[resultId];
        const auto& chunk = batchChunks[resultId];
        std::cerr << "[resultId = " << resultId << "]\n";

        // One chunk can have multiple queries (subreads).
        for (size_t qId = 0; qId < result.size(); ++qId) {
            std::cerr << "    [qId = " << qId << "]\n";

            // Prepare the query data in fwd and rev.
            const char* qSeqFwd = chunk.querySeqs[qId].c_str();
            const int32_t qLen = chunk.querySeqs[qId].size();
            const std::string qSeqRevString = PacBio::Pancake::ReverseComplement(qSeqFwd, 0, qLen);
            const char* qSeqRev = qSeqRevString.c_str();

            // Each query can have multiple mappings.
            for (size_t mapId = 0; mapId < result[qId].mappings.size(); ++mapId) {
                const auto& mapping = result[qId].mappings[mapId];

                std::cerr << "        [mapId = " << mapId << "] "
                          << OverlapWriterBase::PrintOverlapAsM4(mapping->mapping, "", "", true,
                                                                 true)
                          << "\n";

                const char* tSeq = chunk.targetSeqs[mapping->mapping->Bid].c_str();
                const int32_t tLen = chunk.targetSeqs[mapping->mapping->Bid].size();

                // Each mapping is split into regions in between seed hits for alignment.
                std::string qSubSeq;
                std::string tSubSeq;
                for (size_t regId = 0; regId < mapping->regionsForAln.size(); ++regId) {
                    const auto& region = mapping->regionsForAln[regId];
                    std::cerr << "            [regId = " << regId << "] " << region << "\n";

                    if (region.type == RegionType::GLOBAL &&
                            regionsToAdd == BatchAlignerRegionType::SEMIGLOBAL ||
                        region.type != RegionType::GLOBAL &&
                            regionsToAdd == BatchAlignerRegionType::GLOBAL) {
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

std::vector<std::vector<MapperBaseResult>> MapperBatch::StitchAlignments_(
    const std::vector<MapperBatchChunk>& batchChunks,
    const std::vector<std::vector<MapperBaseResult>>& mappingResults,
    const std::vector<AlignmentResult>& internalAlns, const std::vector<AlignmentResult>& flankAlns)
{
    return {};
}

}  // namespace Pancake
}  // namespace PacBio
