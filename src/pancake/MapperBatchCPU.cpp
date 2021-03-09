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

MapperBatchCPU::MapperBatchCPU(const MapperCLRAlignSettings& alignSettings, int32_t numThreads)
    : MapperBatchCPU{alignSettings, nullptr}
{
    fafFallback_ = std::make_unique<Parallel::FireAndForget>(numThreads);
    faf_ = fafFallback_.get();
}

MapperBatchCPU::MapperBatchCPU(const MapperCLRAlignSettings& alignSettings,
                               Parallel::FireAndForget* faf)
    : alignSettings_{alignSettings}, faf_{faf}, fafFallback_{nullptr}
{
}

MapperBatchCPU::~MapperBatchCPU()
{
    if (fafFallback_) {
        fafFallback_->Finalize();
    }
}

std::vector<std::vector<MapperBaseResult>> MapperBatchCPU::MapAndAlign(
    const std::vector<MapperBatchChunk>& batchData)
{
    return MapAndAlignImpl_(batchData, alignSettings_, faf_);
}

void VerboseBatchMappingResults(std::ostream& os,
                                const std::vector<std::vector<MapperBaseResult>>& results)
{
    for (size_t chunkId = 0; chunkId < results.size(); ++chunkId) {
        const auto& chunkResults = results[chunkId];
        for (size_t qId = 0; qId < chunkResults.size(); ++qId) {
            // std::cerr << results[chunkId][qId] << "\n";

            for (size_t mapId = 0; mapId < chunkResults[qId].mappings.size(); ++mapId) {
                if (chunkResults[qId].mappings[mapId] == nullptr) {
                    continue;
                }
                if (chunkResults[qId].mappings[mapId]->mapping == nullptr) {
                    continue;
                }
                const auto& aln = chunkResults[qId].mappings[mapId]->mapping;
                os << *aln << "\n";
            }
        }
    }
}

std::vector<std::vector<MapperBaseResult>> MapperBatchCPU::MapAndAlignImpl_(
    const std::vector<MapperBatchChunk>& batchChunks, const MapperCLRAlignSettings& alignSettings,
    Parallel::FireAndForget* faf)
{
    // Determine how many records should land in each thread, spread roughly evenly.
    const int32_t numRecords = batchChunks.size();
    const int32_t numThreads = faf ? faf->NumThreads() : 1;
    const std::vector<std::pair<int32_t, int32_t>> jobsPerThread =
        PacBio::Pancake::DistributeJobLoad<int32_t>(numThreads, numRecords);

    // Run the mapping in parallel.
    std::vector<std::vector<MapperBaseResult>> results(batchChunks.size());
    const auto Submit = [&jobsPerThread, &batchChunks, &results](int32_t i) {
        const int32_t jobStart = jobsPerThread[i].first;
        const int32_t jobEnd = jobsPerThread[i].second;
        WorkerMapper_(batchChunks, jobStart, jobEnd, results);
    };

    Parallel::Dispatch(faf, jobsPerThread.size(), Submit);

    std::cerr << "[MapperBatchCPU::MapAndAlignImpl_] Mapped:\n";
    VerboseBatchMappingResults(std::cerr, results);

    if (alignSettings.align) {
        // Compute the reverse complements for alignment.
        std::vector<std::vector<std::string>> querySeqsRev =
            ComputeReverseComplements(batchChunks, results, faf);

        // Prepare the sequences for alignment.
        std::vector<PairForBatchAlignment> partsGlobal;
        std::vector<PairForBatchAlignment> partsSemiglobal;
        std::vector<AlignmentStitchInfo> alnStitchInfo;
        int32_t longestSequenceForAln = 0;
        PrepareSequencesForBatchAlignment(batchChunks, querySeqsRev, results,
                                          alignSettings.selfHitPolicy, partsGlobal, partsSemiglobal,
                                          alnStitchInfo, longestSequenceForAln);
        PBLOG_TRACE << "partsGlobal.size() = " << partsGlobal.size();
        PBLOG_TRACE << "partsSemiglobal.size() = " << partsSemiglobal.size();

        // Internal alignment on CPU.
        std::vector<AlignmentResult> internalAlns;
        AlignPartsOnCpu(alignSettings.alignerTypeGlobal, alignSettings.alnParamsGlobal,
                        alignSettings.alignerTypeExt, alignSettings.alnParamsExt, partsGlobal, faf,
                        internalAlns);
        PBLOG_TRACE << "internalAlns.size() = " << internalAlns.size();

        // Flank alignment on CPU.
        std::vector<AlignmentResult> flankAlns;
        AlignPartsOnCpu(alignSettings.alignerTypeGlobal, alignSettings.alnParamsGlobal,
                        alignSettings.alignerTypeExt, alignSettings.alnParamsExt, partsSemiglobal,
                        faf, flankAlns);
        PBLOG_TRACE << "flankAlns.size() = " << flankAlns.size();

        StitchAlignmentsInParallel(results, batchChunks, querySeqsRev, internalAlns, flankAlns,
                                   alnStitchInfo, faf);

        std::cerr << "[MapperBatchCPU::MapAndAlignImpl_] Stitched:\n";
        VerboseBatchMappingResults(std::cerr, results);

        SetUnalignedAndMockedMappings(
            results, alignSettings.selfHitPolicy == MapperSelfHitPolicy::PERFECT_ALIGNMENT,
            alignSettings.alnParamsGlobal.matchScore);

        std::cerr << "[MapperBatchCPU::MapAndAlignImpl_] SetUnalignedAndMockedMappings:\n";
        VerboseBatchMappingResults(std::cerr, results);

        UpdateSecondaryAndFilter(results, faf, batchChunks);

        std::cerr << "[MapperBatchCPU::MapAndAlignImpl_] UpdateSecondaryAndFilter:\n";
        VerboseBatchMappingResults(std::cerr, results);
    }

    return results;
}

void MapperBatchCPU::WorkerMapper_(const std::vector<MapperBatchChunk>& batchChunks,
                                   int32_t startId, int32_t endId,
                                   std::vector<std::vector<MapperBaseResult>>& results)
{
    for (int32_t i = startId; i < endId; ++i) {
        const auto& chunk = batchChunks[i];

        // Create a copy of the settings so that we can turn off the alignment.
        MapperCLRSettings settingsCopy;
        settingsCopy.map = chunk.mapSettings;
        settingsCopy.align.align = false;

        // Create the mapper.
        MapperCLR mapper(settingsCopy);

        results[i] = mapper.MapAndAlign(chunk.targetSeqs, chunk.querySeqs);
    }
}

void UpdateSecondaryAndFilter(std::vector<std::vector<MapperBaseResult>>& mappingResults,
                              Parallel::FireAndForget* faf,
                              const std::vector<MapperBatchChunk>& batchChunks)
{

    const int32_t numThreads = faf ? faf->NumThreads() : 1;
    const int32_t numEntries = mappingResults.size();
    const std::vector<std::pair<int32_t, int32_t>> jobsPerThread =
        PacBio::Pancake::DistributeJobLoad<int32_t>(numThreads, numEntries);

    // Results are a vector for every chunk (one chunk is one ZMW).
    const auto Submit = [&](int32_t jobId) {
        std::cerr << "jobId = " << jobId << " / " << jobsPerThread.size() << "\n";
        const int32_t jobStart = jobsPerThread[jobId].first;
        const int32_t jobEnd = jobsPerThread[jobId].second;

        std::cerr << "jobStart = " << jobStart << ", jobEnd = " << jobEnd << "\n";
        std::cerr << "\n";

        for (int32_t resultId = jobStart; resultId < jobEnd; ++resultId) {
            const auto& settings = batchChunks[resultId].mapSettings;
            auto& result = mappingResults[resultId];
            // One chunk can have multiple queries (subreads).
            for (size_t qId = 0; qId < result.size(); ++qId) {
                std::cerr << "[" << __FUNCTION__ << "] Before WrapFlagSecondaryAndSupplementary:\n";
                std::cerr << result[qId] << "\n";
                // for (size_t mapId = 0; mapId < result[qId].mappings.size(); ++mapId) {
                //     if (result[qId].mappings[mapId] == nullptr) {
                //         std::cerr << "[resultId = " << resultId << ", qId = " << qId << ", mapId = " << mapId << "] result[qId].mappings[mapId] == nullptr\n";
                //         continue;
                //     }
                //     if (result[qId].mappings[mapId]->mapping == nullptr) {
                //         std::cerr << "[resultId = " << resultId << ", qId = " << qId << ", mapId = " << mapId << "] result[qId].mappings[mapId]->mapping == nullptr\n";
                //         continue;
                //     }
                //     const auto& aln = result[qId].mappings[mapId]->mapping;
                //     std::cerr << "[resultId = " << resultId << ", qId = " << qId << ", mapId = " << mapId << "] " << *aln << "\n";
                // }

                // Secondary/supplementary flagging.
                WrapFlagSecondaryAndSupplementary(result[qId].mappings,
                                                  settings.secondaryAllowedOverlapFractionQuery,
                                                  settings.secondaryAllowedOverlapFractionTarget,
                                                  settings.secondaryMinScoreFraction);

                std::cerr << "[" << __FUNCTION__ << "] After WrapFlagSecondaryAndSupplementary:\n";
                std::cerr << result[qId] << "\n";
                // for (size_t mapId = 0; mapId < result[qId].mappings.size(); ++mapId) {
                //     if (result[qId].mappings[mapId] == nullptr) {
                //         std::cerr << "[resultId = " << resultId << ", qId = " << qId << ", mapId = " << mapId << "] result[qId].mappings[mapId] == nullptr\n";
                //         continue;
                //     }
                //     if (result[qId].mappings[mapId]->mapping == nullptr) {
                //         std::cerr << "[resultId = " << resultId << ", qId = " << qId << ", mapId = " << mapId << "] result[qId].mappings[mapId]->mapping == nullptr\n";
                //         continue;
                //     }
                //     c