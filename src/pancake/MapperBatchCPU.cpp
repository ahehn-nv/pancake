// Authors: Ivan Sovic

#include <pacbio/alignment/AlignmentTools.h>
#include <pacbio/pancake/AlignerBase.h>
#include <pacbio/pancake/MapperBatchCPU.h>
#include <pacbio/pancake/MapperBatchUtility.h>
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

MapperBatchCPU::MapperBatchCPU(const MapperCLRSettings& settings, int32_t numThreads)
    : MapperBatchCPU{settings, nullptr}
{
    fafFallback_ = std::make_unique<Parallel::FireAndForget>(numThreads);
    faf_ = fafFallback_.get();
}

MapperBatchCPU::MapperBatchCPU(const MapperCLRSettings& settings, Parallel::FireAndForget* faf)
    : settings_{settings}, faf_{faf}, fafFallback_{nullptr}
{
}

MapperBatchCPU::~MapperBatchCPU()
{
    if (fafFallback_) {
        fafFallback_->Finalize();
    }
}

std::vector<std::vector<MapperBaseResult>> MapperBatchCPU::DummyMapAndAlign(
    const std::vector<MapperBatchChunk>& batchData)
{
    return DummyMapAndAlignImpl_(batchData, settings_, faf_);
}

std::vector<std::vector<MapperBaseResult>> MapperBatchCPU::MapAndAlign(
    const std::vector<MapperBatchChunk>& batchData)
{
    return MapAndAlignImpl_(batchData, settings_, faf_);
}

std::vector<std::vector<MapperBaseResult>> MapperBatchCPU::DummyMapAndAlignImpl_(
    const std::vector<MapperBatchChunk>& batchChunks, MapperCLRSettings settings,
    Parallel::FireAndForget* /*faf*/)
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

std::vector<std::vector<MapperBaseResult>> MapperBatchCPU::MapAndAlignImpl_(
    const std::vector<MapperBatchChunk>& batchChunks, MapperCLRSettings settings,
    Parallel::FireAndForget* faf)
{
    // Determine how many records should land in each thread, spread roughly evenly.
    const int32_t numRecords = batchChunks.size();
    const int32_t numThreads = faf ? faf->NumThreads() : 1;
    const std::vector<std::pair<int32_t, int32_t>> jobsPerThread =
        PacBio::Pancake::DistributeJobLoad<int32_t>(numThreads, numRecords);

    // Create a mapper for each thread.
    MapperCLRSettings settingsCopy = settings;
    settingsCopy.align = false;
    std::vector<std::unique_ptr<MapperCLR>> mappers;
    for (size_t i = 0; i < jobsPerThread.size(); ++i) {
        auto mapper = std::make_unique<MapperCLR>(settingsCopy);
        mappers.emplace_back(std::move(mapper));
    }

    // Run the mapping in parallel.
    std::vector<std::vector<MapperBaseResult>> results(batchChunks.size());
    const auto Submit = [&jobsPerThread, &batchChunks, &mappers, &results](int32_t i) {
        const int32_t jobStart = jobsPerThread[i].first;
        const int32_t jobEnd = jobsPerThread[i].second;
        WorkerMapper_(batchChunks, jobStart, jobEnd, mappers[i], results);
    };

    Parallel::Dispatch(faf, jobsPerThread.size(), Submit);

    if (settings.align) {
        // Compute the reverse complements for alignment.
        std::vector<std::vector<std::string>> querySeqsRev;
        for (const auto& chunk : batchChunks) {
            std::vector<std::string> revSeqs;
            for (const auto& query : chunk.querySeqs) {
                revSeqs.emplace_back(
                    PacBio::Pancake::ReverseComplement(query.c_str(), 0, query.size()));
            }
            querySeqsRev.emplace_back(std::move(revSeqs));
        }

        // Prepare the sequences for alignment.
        std::vector<PairForBatchAlignment> partsGlobal;
        std::vector<PairForBatchAlignment> partsSemiglobal;
        std::vector<AlignmentStitchInfo> alnStitchInfo;
        int32_t longestSequenceForAln = 0;
        PrepareSequencesForBatchAlignment(batchChunks, querySeqsRev, results, partsGlobal,
                                          partsSemiglobal, alnStitchInfo, longestSequenceForAln);
        PBLOG_TRACE << "partsGlobal.size() = " << partsGlobal.size();
        PBLOG_TRACE << "partsSemiglobal.size() = " << partsSemiglobal.size();

        // Internal alignment on CPU.
        std::vector<AlignmentResult> internalAlns;
        AlignPartsOnCpu(settings.alignerTypeGlobal, settings.alnParamsGlobal,
                        settings.alignerTypeExt, settings.alnParamsExt, partsGlobal, faf,
                        internalAlns);
        PBLOG_TRACE << "internalAlns.size() = " << internalAlns.size();

        // Flank alignment on CPU.
        std::vector<AlignmentResult> flankAlns;
        AlignPartsOnCpu(settings.alignerTypeGlobal, settings.alnParamsGlobal,
                        settings.alignerTypeExt, settings.alnParamsExt, partsSemiglobal, faf,
                        flankAlns);
        PBLOG_TRACE << "flankAlns.size() = " << flankAlns.size();

        StitchAlignments(results, batchChunks, querySeqsRev, internalAlns, flankAlns,
                         alnStitchInfo);

        UpdateSecondaryAndFilter(results, settings.secondaryAllowedOverlapFractionQuery,
                                 settings.secondaryAllowedOverlapFractionTarget,
                                 settings.secondaryMinScoreFraction, settings.bestNSecondary);
    }

    return results;
}

void MapperBatchCPU::WorkerMapper_(const std::vector<MapperBatchChunk>& batchChunks,
                                   int32_t startId, int32_t endId,
                                   const std::unique_ptr<MapperCLR>& mapper,
                                   std::vector<std::vector<MapperBaseResult>>& results)
{
    for (int32_t i = startId; i < endId; ++i) {
        const auto& chunk = batchChunks[i];
        results[i] = mapper->MapAndAlign(chunk.targetSeqs, chunk.querySeqs);
    }
}

void UpdateSecondaryAndFilter(std::vector<std::vector<MapperBaseResult>>& mappingResults,
                              double secondaryAllowedOverlapFractionQuery,
                              double secondaryAllowedOverlapFractionTarget,
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

int32_t AlignPartsOnCpu(const AlignerType& alignerTypeGlobal,
                        const AlignmentParameters& alnParamsGlobal,
                        const AlignerType& alignerTypeExt, const AlignmentParameters& alnParamsExt,
                        const std::vector<PairForBatchAlignment>& parts, const int32_t numThreads,
                        std::vector<AlignmentResult>& retAlns)
{
    Parallel::FireAndForget faf(numThreads);
    const int32_t result = AlignPartsOnCpu(alignerTypeGlobal, alnParamsGlobal, alignerTypeExt,
                                           alnParamsExt, parts, &faf, retAlns);
    faf.Finalize();
    return result;
}
int32_t AlignPartsOnCpu(const AlignerType& alignerTypeGlobal,
                        const AlignmentParameters& alnParamsGlobal,
                        const AlignerType& alignerTypeExt, const AlignmentParameters& alnParamsExt,
                        const std::vector<PairForBatchAlignment>& parts,
                        Parallel::FireAndForget* faf, std::vector<AlignmentResult>& retAlns)
{
    retAlns.resize(parts.size());

    std::vector<size_t> partIds;

    AlignerBatchCPU aligner(alignerTypeGlobal, alnParamsGlobal, alignerTypeExt, alnParamsExt);

    for (size_t i = 0; i < parts.size(); ++i) {
        const auto& part = parts[i];
        if (retAlns[i].valid) {
            continue;
        }
        partIds.emplace_back(i);
        if (part.regionType == RegionType::FRONT) {
            // Reverse the sequences for front flank alignment. No need to complement.
            std::string query(part.query, part.queryLen);
            std::reverse(query.begin(), query.end());
            std::string target(part.target, part.targetLen);
            std::reverse(target.begin(), target.end());
            aligner.AddSequencePair(query.c_str(), part.queryLen, target.c_str(), part.targetLen,
                                    part.regionType == RegionType::GLOBAL);
        } else {
            aligner.AddSequencePair(part.query, part.queryLen, part.target, part.targetLen,
                                    part.regionType == RegionType::GLOBAL);
        }
    }
    aligner.AlignAll(faf);

    const std::vector<AlignmentResult>& partInternalAlns = aligner.GetAlnResults();
    int32_t numNotValid = 0;
    for (size_t i = 0; i < partInternalAlns.size(); ++i) {
        const auto& aln = partInternalAlns[i];
        if (aln.valid == false) {
            ++numNotValid;
        }
        retAlns[partIds[i]] = std::move(partInternalAlns[i]);
    }
    return numNotValid;
}

}  // namespace Pancake
}  // namespace PacBio
