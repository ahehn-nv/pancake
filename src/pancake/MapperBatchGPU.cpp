// Authors: Ivan Sovic
#include <pacbio/pancake/AlignerBatchGPU.h>
#include <pacbio/pancake/MapperBatchGPU.h>

#include <pacbio/alignment/AlignmentTools.h>
#include <pacbio/pancake/AlignerBase.h>
#include <pacbio/pancake/AlignmentSeeded.h>
#include <pacbio/pancake/OverlapWriterFactory.h>
#include <pacbio/util/TicToc.h>
#include <pbcopper/logging/Logging.h>
#include <pbcopper/parallel/FireAndForget.h>
#include <pbcopper/parallel/WorkQueue.h>
#include <pbcopper/utility/Stopwatch.h>
#include <algorithm>
#include <array>
#include <cassert>
#include <claraparabricks/genomeworks/utils/allocator.hpp>
#include <iostream>
#include <sstream>
#include <tuple>

namespace PacBio {
namespace Pancake {

MapperBatchGPU::MapperBatchGPU(const MapperCLRSettings& settings, int32_t numThreads,
                               int32_t gpuStartBandwidth, int32_t gpuMaxBandwidth,
                               uint32_t gpuDeviceId, int64_t gpuMemoryBytes,
                               bool alignRemainingOnCpu)
    : settings_{settings}
    , gpuStartBandwidth_(gpuStartBandwidth)
    , gpuMaxBandwidth_(gpuMaxBandwidth)
    , alignRemainingOnCpu_(alignRemainingOnCpu)
    , faf_{nullptr}
    , fafFallback_(nullptr)
    , aligner_{nullptr}
{
    fafFallback_ = std::make_unique<Parallel::FireAndForget>(numThreads);
    faf_ = fafFallback_.get();
    aligner_ = std::make_unique<AlignerBatchGPU>(faf_, settings.alnParamsGlobal, gpuStartBandwidth,
                                                 gpuDeviceId, gpuMemoryBytes);
}

MapperBatchGPU::MapperBatchGPU(const MapperCLRSettings& settings, Parallel::FireAndForget* faf,
                               int32_t gpuStartBandwidth, int32_t gpuMaxBandwidth,
                               uint32_t gpuDeviceId, int64_t gpuMemoryBytes,
                               bool alignRemainingOnCpu)
    : settings_{settings}
    , gpuStartBandwidth_(gpuStartBandwidth)
    , gpuMaxBandwidth_(gpuMaxBandwidth)
    , alignRemainingOnCpu_(alignRemainingOnCpu)
    , faf_{faf}
    , fafFallback_(nullptr)
    , aligner_{std::make_unique<AlignerBatchGPU>(faf, settings.alnParamsGlobal, gpuStartBandwidth,
                                                 gpuDeviceId, gpuMemoryBytes)}
{}

MapperBatchGPU::~MapperBatchGPU()
{
    if (fafFallback_) {
        fafFallback_->Finalize();
    }
}

std::vector<std::vector<MapperBaseResult>> MapperBatchGPU::MapAndAlign(
    const std::vector<MapperBatchChunk>& batchData)
{
    assert(aligner_);
    return MapAndAlignImpl_(batchData, settings_, alignRemainingOnCpu_, gpuStartBandwidth_,
                            gpuMaxBandwidth_, *aligner_, faf_);
}

std::vector<std::vector<std::string>> ComputeReverseComplements(
    const std::vector<MapperBatchChunk>& batchChunks,
    const std::vector<std::vector<MapperBaseResult>>& mappingResults, Parallel::FireAndForget* faf)
{
    PacBio::Utility::Stopwatch timer;

    std::vector<std::vector<uint8_t>> shouldReverse(batchChunks.size());
    for (size_t i = 0; i < batchChunks.size(); ++i) {
        shouldReverse[i].resize(batchChunks[i].querySeqs.size(), false);
    }

    // Results are a vector for every chunk (one chunk is one ZMW).
    for (size_t chunkId = 0; chunkId < mappingResults.size(); ++chunkId) {
        auto& result = mappingResults[chunkId];
        // One chunk can have multiple queries (subreads).
        for (size_t qId = 0; qId < result.size(); ++qId) {
            for (size_t mapId = 0; mapId < result[qId].mappings.size(); ++mapId) {
                if (result[qId].mappings[mapId] == nullptr ||
                    result[qId].mappings[mapId]->mapping == nullptr) {
                    continue;
                }
                const OverlapPtr& aln = result[qId].mappings[mapId]->mapping;
                shouldReverse[chunkId][qId] |= aln->Brev;
            }
        }
    }

    // Determine how many records should land in each thread, spread roughly evenly.
    const int32_t numThreads = faf ? faf->NumThreads() : 1;
    const int32_t numRecords = batchChunks.size();
    const std::vector<std::pair<int32_t, int32_t>> jobsPerThread =
        PacBio::Pancake::DistributeJobLoad<int32_t>(numThreads, numRecords);

    std::vector<std::vector<std::string>> querySeqsRev(batchChunks.size());

    timer.Freeze();
    PBLOG_INFO << "CPU RevComp prepare    : " << timer.ElapsedTime();
    timer.Reset();

    const auto Submit = [&batchChunks, &jobsPerThread, &shouldReverse, &querySeqsRev](int32_t idx) {
        const int32_t jobStart = jobsPerThread[idx].first;
        const int32_t jobEnd = jobsPerThread[idx].second;
        for (int32_t chunkId = jobStart; chunkId < jobEnd; ++chunkId) {
            auto& revSeqs = querySeqsRev[chunkId];
            for (size_t qId = 0; qId < batchChunks[chunkId].querySeqs.size(); ++qId) {
                if (shouldReverse[chunkId][qId]) {
                    const auto& query = batchChunks[chunkId].querySeqs[qId];
                    revSeqs.emplace_back(
                        PacBio::Pancake::ReverseComplement(query.c_str(), 0, query.size()));
                } else {
                    revSeqs.emplace_back("");
                }
            }
        }
    };
    Parallel::Dispatch(faf, jobsPerThread.size(), Submit);

    timer.Freeze();
    PBLOG_INFO << "CPU RevComp compute    : " << timer.ElapsedTime();
    timer.Reset();

    return querySeqsRev;
}

std::vector<std::vector<MapperBaseResult>> MapperBatchGPU::MapAndAlignImpl_(
    const std::vector<MapperBatchChunk>& batchChunks, const MapperCLRSettings& settings,
    bool alignRemainingOnCpu, int32_t gpuStartBandwidth, int32_t gpuMaxBandwidth,
    AlignerBatchGPU& aligner, Parallel::FireAndForget* faf)
{
    PacBio::Utility::Stopwatch timer;
    int64_t cpuTime = 0;
    int64_t gpuTime = 0;

    const int32_t numThreads = faf ? faf->NumThreads() : 1;
    // Determine how many records should land in each thread, spread roughly evenly.
    const int32_t numRecords = batchChunks.size();
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

    const auto Submit = [&batchChunks, &mappers, &jobsPerThread, &results](int32_t idx) {
        const int32_t jobStart = jobsPerThread[idx].first;
        const int32_t jobEnd = jobsPerThread[idx].second;
        WorkerMapper_(batchChunks, jobStart, jobEnd, *mappers[idx], results);
    };
    Parallel::Dispatch(faf, jobsPerThread.size(), Submit);
    timer.Freeze();
    cpuTime += timer.ElapsedNanoseconds();
    PBLOG_INFO << "CPU Mapping            : " << timer.ElapsedTime();

    // Align the mappings if required.
    if (settings.align) {
        int64_t seedAlnCpuTime = 0;
        int64_t seedAlnGpuTime = 0;
        timer.Reset();
        // Sanity check.
        if (gpuStartBandwidth <= 0) {
            throw std::runtime_error(
                "The startBandwidth needs to be a positive number. gpuStartBandwidth = " +
                std::to_string(gpuStartBandwidth));
        }

        // Compute the reverse complements for alignment.
        std::vector<std::vector<std::string>> querySeqsRev =
            ComputeReverseComplements(batchChunks, results, faf);
        timer.Freeze();
        seedAlnCpuTime += timer.ElapsedNanoseconds();
        PBLOG_INFO << "CPU RevComp            : " << timer.ElapsedTime();
        timer.Reset();

        PBLOG_TRACE << "Preparing parts for alignment.";

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

        // Global alignment on GPU. Try using different bandwidths,
        // increasing the bandwidth for the failed parts each iteration.
        std::vector<AlignmentResult> internalAlns;
        int32_t currentBandwidth = gpuStartBandwidth;
        int32_t numInternalNotValid = 0;
        const int32_t maxBandwidth =
            (gpuMaxBandwidth <= 0) ? (longestSequenceForAln * 2 + 1) : gpuMaxBandwidth;

        while (true) {
            PBLOG_TRACE << "Trying bandwidth: " << currentBandwidth;
            timer.Reset();
            aligner.ResetMaxBandwidth(gpuMaxBandwidth);
            timer.Freeze();
            cpuTime += timer.ElapsedNanoseconds();
            PBLOG_INFO << "CPU Reset BW           : " << timer.ElapsedTime();
            timer.Reset();
            numInternalNotValid = AlignPartsOnGPU_(aligner, partsGlobal, internalAlns,
                                                   seedAlnCpuTime, seedAlnGpuTime);
            if (numInternalNotValid == 0) {
                break;
            }
            if (currentBandwidth >= maxBandwidth) {
                break;
            }
            currentBandwidth *= 2;
            // Ensure that the last bandwidth tried is the maxBandwidth.
            if (currentBandwidth > maxBandwidth) {
                currentBandwidth = maxBandwidth;
            }
        }
        cpuTime += seedAlnCpuTime;
        gpuTime += seedAlnGpuTime;
        PBLOG_INFO << "CPU Internal Prepare   : "
                   << PacBio::Utility::Stopwatch::PrettyPrintNanoseconds(seedAlnCpuTime);
        PBLOG_INFO << "GPU Internal Alignment : "
                   << PacBio::Utility::Stopwatch::PrettyPrintNanoseconds(seedAlnGpuTime);
        timer.Reset();

        // Fallback to the CPU if there are any unaligned parts left.
        int64_t prepareTime = 0;
        int64_t alignTime = 0;
        if (alignRemainingOnCpu && numInternalNotValid > 0) {
            PBLOG_TRACE << "Trying to align remaining parts on CPU.";
            const int32_t numNotValidInternal = AlignPartsOnCpu(
                settings.alignerTypeGlobal, settings.alnParamsGlobal, settings.alignerTypeExt,
                settings.alnParamsExt, partsGlobal, faf, internalAlns, prepareTime, alignTime);
            PBLOG_TRACE << "Total not valid: " << numNotValidInternal << " / "
                        << internalAlns.size() << "\n";
        }
        PBLOG_TRACE << "internalAlns.size() = " << internalAlns.size();
        timer.Freeze();
        cpuTime += timer.ElapsedNanoseconds();
        PBLOG_INFO << "CPU Fallback           : " << timer.ElapsedTime() << ' '
                   << numInternalNotValid;
        timer.Reset();

        // Flank alignment on CPU.
        std::vector<AlignmentResult> flankAlns;
        prepareTime = 0;
        alignTime = 0;
        const int32_t numNotValidFlanks = AlignPartsOnCpu(
            settings.alignerTypeGlobal, settings.alnParamsGlobal, settings.alignerTypeExt,
            settings.alnParamsExt, partsSemiglobal, faf, flankAlns, prepareTime, alignTime);
        PBLOG_TRACE << "Total not valid: " << numNotValidFlanks << " / " << flankAlns.size()
                    << "\n";
        timer.Freeze();
        cpuTime += timer.ElapsedNanoseconds();
        PBLOG_INFO << "CPU Flanks Prepare     : "
                   << PacBio::Utility::Stopwatch::PrettyPrintNanoseconds(prepareTime);
        PBLOG_INFO << "CPU Flanks Alignment   : "
                   << PacBio::Utility::Stopwatch::PrettyPrintNanoseconds(alignTime);
        timer.Reset();

        StitchAlignmentsInParallel(results, batchChunks, querySeqsRev, internalAlns, flankAlns,
                                   alnStitchInfo, faf);
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
    PBLOG_INFO << "GPU Time               : "
               << PacBio::Utility::Stopwatch::PrettyPrintNanoseconds(gpuTime);

    return results;
}

void MapperBatchGPU::WorkerMapper_(const std::vector<MapperBatchChunk>& batchChunks,
                                   int32_t startId, int32_t endId, MapperCLR& mapper,
                                   std::vector<std::vector<MapperBaseResult>>& results)
{
    for (int32_t i = startId; i < endId; ++i) {
        const auto& chunk = batchChunks[i];
        results[i] = mapper.MapAndAlign(chunk.targetSeqs, chunk.querySeqs);
    }
}

int32_t MapperBatchGPU::AlignPartsOnGPU_(AlignerBatchGPU& aligner,
                                         const std::vector<PairForBatchAlignment>& parts,
                                         std::vector<AlignmentResult>& retInternalAlns,
                                         int64_t& cpuTime, int64_t& gpuTime)
{
    PacBio::Utility::Stopwatch timer;
    retInternalAlns.resize(parts.size());
    cpuTime += timer.ElapsedNanoseconds();

    int32_t totalNumNotValid = 0;
    size_t partId = 0;
    while (partId < parts.size()) {
        timer.Reset();
        aligner.Clear();

        std::vector<size_t> partIds;

        PBLOG_TRACE << "Preparing sequences for GPU alignment.";
        for (; partId < parts.size(); ++partId) {
            const PairForBatchAlignment& part = parts[partId];

            if (retInternalAlns[partId].valid) {
                continue;
            }
            partIds.emplace_back(partId);

            StatusAddSequencePair rv;
            if (part.regionType == RegionType::FRONT) {
                // Reverse the sequences for front flank alignment. No need to complement.
                std::string query(part.query, part.queryLen);
                std::reverse(query.begin(), query.end());
                std::string target(part.target, part.targetLen);
                std::reverse(target.begin(), target.end());
                rv = aligner.AddSequencePair(query.c_str(), part.queryLen, target.c_str(),
                                             part.targetLen);
            } else {
                rv =
                    aligner.AddSequencePair(part.query, part.queryLen, part.target, part.targetLen);
            }

            if (rv == StatusAddSequencePair::EXCEEDED_MAX_ALIGNMENTS) {
                break;

            } else if (rv != StatusAddSequencePair::OK) {
                throw std::runtime_error(
                    "Error occurred while trying to add sequences for batch alignment to "
                    "AlignerBatchGPU.");
            }
        }
        cpuTime += timer.ElapsedNanoseconds();

        PBLOG_TRACE << "Aligning batch of " << aligner.BatchSize() << " sequence pairs.";
        std::pair<int64_t, int64_t> seedTimings = aligner.AlignAll();
        cpuTime += seedTimings.first;
        gpuTime += seedTimings.second;
        timer.Reset();

        const std::vector<AlignmentResult>& partInternalAlns = aligner.GetAlnResults();

        int32_t numNotValid = 0;
        for (size_t i = 0; i < partInternalAlns.size(); ++i) {
            const auto& aln = partInternalAlns[i];
            if (aln.valid == false) {
                ++numNotValid;
            }
            // retInternalAlns.emplace_back(std::move(partInternalAlns[i]));
            retInternalAlns[partIds[i]] = std::move(partInternalAlns[i]);
        }
        totalNumNotValid += numNotValid;
        cpuTime += timer.ElapsedNanoseconds();
    }
    PBLOG_TRACE << "Total not valid: " << totalNumNotValid << " / " << retInternalAlns.size()
                << "\n";
    return totalNumNotValid;
}

}  // namespace Pancake
}  // namespace PacBio
