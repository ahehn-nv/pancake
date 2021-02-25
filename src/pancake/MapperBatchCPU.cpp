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
#include <pbcopper/utility/Stopwatch.h>
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
{}

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
    PacBio::Utility::Stopwatch timer;
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
    PBLOG_INFO << "CPU Mapping            : " << timer.ElapsedTime();
    timer.Reset();

    int64_t cpuTime = 0;
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
        timer.Freeze();
        cpuTime += timer.ElapsedNanoseconds();
        PBLOG_INFO << "CPU RevComp            : " << timer.ElapsedTime();
        timer.Reset();

        // Prepare the sequences for alignment.
        std::vector<PairForBatchAlignment> partsGlobal;
        std::vector<PairForBatchAlignment> partsSemiglobal;
        std::vector<AlignmentStitchInfo> alnStitchInfo;
        int32_t longestSequenceForAln = 0;
        PrepareSequencesForBatchAlignment(batchChunks, querySeqsRev, results, partsGlobal,
                                          partsSemiglobal, alnStitchInfo, longestSequenceForAln);
        PBLOG_TRACE << "partsGlobal.size() = " << partsGlobal.size();
        PBLOG_TRACE << "partsSemiglobal.size() = " << partsSemiglobal.size();
        timer.Freeze();
        cpuTime += timer.ElapsedNanoseconds();
        PBLOG_INFO << "CPU Prepare            : " << timer.ElapsedTime();
        timer.Reset();

        // Internal alignment on CPU.
        std::vector<AlignmentResult> internalAlns;
        int64_t prepareTime = 0;
        int64_t alignTime = 0;
        AlignPartsOnCpu(settings.alignerTypeGlobal, settings.alnParamsGlobal,
                        settings.alignerTypeExt, settings.alnParamsExt, partsGlobal, faf,
                        internalAlns, prepareTime, alignTime);
        PBLOG_TRACE << "internalAlns.size() = " << internalAlns.size();
        timer.Freeze();
        cpuTime += timer.ElapsedNanoseconds();
        PBLOG_INFO << "CPU Internal Prepare   : "
                   << PacBio::Utility::Stopwatch::PrettyPrintNanoseconds(prepareTime);
        PBLOG_INFO << "CPU Internal Alignment : "
                   << PacBio::Utility::Stopwatch::PrettyPrintNanoseconds(alignTime);
        timer.Reset();

        // Flank alignment on CPU.
        std::vector<AlignmentResult> flankAlns;
        prepareTime = 0;
        alignTime = 0;
        AlignPartsOnCpu(settings.alignerTypeGlobal, settings.alnParamsGlobal,
                        settings.alignerTypeExt, settings.alnParamsExt, partsSemiglobal, faf,
                        flankAlns, prepareTime, alignTime);
        PBLOG_TRACE << "flankAlns.size() = " << flankAlns.size();
        timer.Freeze();
        cpuTime += timer.ElapsedNanoseconds();
        PBLOG_INFO << "CPU Flanks Prepare     : "
                   << PacBio::Utility::Stopwatch::PrettyPrintNanoseconds(prepareTime);
        PBLOG_INFO << "CPU Flanks Alignment   : "
                   << PacBio::Utility::Stopwatch::PrettyPrintNanoseconds(alignTime);
        timer.Reset();

        StitchAlignments(results, batchChunks, querySeqsRev, internalAlns, flankAlns,
                         alnStitchInfo);
        timer.Freeze();
        cpuTime += timer.ElapsedNanoseconds();
        PBLOG_INFO << "CPU Stitch             : " << timer.ElapsedTime();
        timer.Reset();

        UpdateSecondaryAndFilter(results, settings.secondaryAllowedOverlapFractionQuery,
                                 settings.secondaryAllowedOverlapFractionTarget,
                                 settings.secondaryMinScoreFraction, settings.bestNSecondary, faf);
        timer.Freeze();
        cpuTime += timer.ElapsedNanoseconds();
        PBLOG_INFO << "CPU Update             : " << timer.ElapsedTime();
        timer.Reset();
    }
    PBLOG_INFO << "CPU Time               : "
               << PacBio::Utility::Stopwatch::PrettyPrintNanoseconds(cpuTime);

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
                              double secondaryMinScoreFraction, int32_t bestNSecondary,
                              Parallel::FireAndForget* faf)
{
    // Results are a vector for every chunk (one chunk is one ZMW).
    const auto Submit = [&](int32_t resultId) {
        auto& result = mappingResults[resultId];
        // One chunk can have multiple queries (subreads).
        for (size_t qId = 0; qId < result.size(); ++qId) {
            // Secondary/supplementary flagging.
            WrapFlagSecondaryAndSupplementary(
                result[qId].mappings, secondaryAllowedOverlapFractionQuery,
                secondaryAllowedOverlapFractionTarget, secondaryMinScoreFraction);
            CondenseMappings(result[qId].mappings, bestNSecondary);
        }
    };
    const int32_t numThreads = faf ? faf->NumThreads() : 1;
    const int32_t numEntries = mappingResults.size();
    const std::vector<std::pair<int32_t, int32_t>> jobsPerThread =
        PacBio::Pancake::DistributeJobLoad<int32_t>(numThreads, numEntries);
    Parallel::Dispatch(faf, numEntries, Submit);
}

int32_t AlignPartsOnCpu(const AlignerType& alignerTypeGlobal,
                        const AlignmentParameters& alnParamsGlobal,
                        const AlignerType& alignerTypeExt, const AlignmentParameters& alnParamsExt,
                        const std::vector<PairForBatchAlignment>& parts, const int32_t numThreads,
                        std::vector<AlignmentResult>& retAlns, int64_t& prepareTime,
                        int64_t& alignTime)
{
    Parallel::FireAndForget faf(numThreads);
    const int32_t result =
        AlignPartsOnCpu(alignerTypeGlobal, alnParamsGlobal, alignerTypeExt, alnParamsExt, parts,
                        &faf, retAlns, prepareTime, alignTime);
    faf.Finalize();
    return result;
}
int32_t AlignPartsOnCpu(const AlignerType& alignerTypeGlobal,
                        const AlignmentParameters& alnParamsGlobal,
                        const AlignerType& alignerTypeExt, const AlignmentParameters& alnParamsExt,
                        const std::vector<PairForBatchAlignment>& parts,
                        Parallel::FireAndForget* faf, std::vector<AlignmentResult>& retAlns,
                        int64_t& prepareTime, int64_t& alignTime)
{
    PacBio::Utility::Stopwatch timer;
    retAlns.resize(parts.size());

    std::vector<size_t> partIds;

    AlignerBatchCPU aligner(faf, alignerTypeGlobal, alnParamsGlobal, alignerTypeExt, alnParamsExt);

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
    prepareTime += timer.ElapsedNanoseconds();
    timer.Reset();

    aligner.AlignAll();
    alignTime += timer.ElapsedNanoseconds();
    timer.Reset();

    const std::vector<AlignmentResult>& partInternalAlns = aligner.GetAlnResults();
    int32_t numNotValid = 0;
    for (size_t i = 0; i < partInternalAlns.size(); ++i) {
        const auto& aln = partInternalAlns[i];
        if (aln.valid == false) {
            ++numNotValid;
        }
        retAlns[partIds[i]] = std::move(partInternalAlns[i]);
    }
    prepareTime += timer.ElapsedNanoseconds();
    return numNotValid;
}

}  // namespace Pancake
}  // namespace PacBio
