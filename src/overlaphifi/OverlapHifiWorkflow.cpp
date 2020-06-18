// Authors: Ivan Sovic

#include "overlaphifi/OverlapHifiWorkflow.h"
#include <pacbio/overlaphifi/Mapper.h>
#include <pacbio/overlaphifi/OverlapHifiSettings.h>
#include <pacbio/overlaphifi/OverlapWriterFactory.h>
#include <pacbio/overlaphifi/SeedIndex.h>
#include <pacbio/seeddb/Seed.h>
#include <pacbio/seeddb/SeedDBIndexCache.h>
#include <pacbio/seeddb/SeedDBReaderCachedBlock.h>
#include <pacbio/seeddb/SeedDBReaderRawBlock.h>
#include <pacbio/seqdb/SeqDBIndexCache.h>
#include <pacbio/seqdb/SeqDBReaderCached.h>
#include <pacbio/util/TicToc.h>
#include <pbcopper/logging/LogLevel.h>
#include <pbcopper/logging/Logging.h>
#include <pbcopper/parallel/FireAndForget.h>
#include <pbcopper/parallel/WorkQueue.h>
#include <sstream>

namespace PacBio {
namespace Pancake {

void Worker(const PacBio::Pancake::SeqDBReaderCachedBlock& targetSeqDBReader,
            const PacBio::Pancake::SeedIndex& index,
            const PacBio::Pancake::SeqDBReaderCachedBlock& querySeqDBReader,
            const PacBio::Pancake::SeedDBReaderCachedBlock& querySeedDBReader,
            const OverlapHifiSettings& /*settings*/, const PacBio::Pancake::Mapper& mapper,
            int64_t freqCutoff, bool generateFlippedOverlaps, int32_t start, int32_t end,
            std::vector<PacBio::Pancake::MapperResult>& results)
{
    int32_t numRecords = querySeqDBReader.records().size();
    if (start < 0 || end < 0 || start > end || start > numRecords || end > numRecords) {
        std::ostringstream oss;
        oss << "Invalid start/end indexes provided to the Worker. start = " << start
            << ", end = " << end << ", numRecords = " << numRecords;
        throw std::runtime_error(oss.str());
    }

    // Map the reads.
    for (int32_t i = start; i < end; ++i) {
        const auto& querySeq = querySeqDBReader.records()[i];
        const auto& querySeeds = querySeedDBReader.GetSeedsForSequence(querySeq.Id());
        results[i] = mapper.Map(targetSeqDBReader, index, querySeq, querySeeds, freqCutoff,
                                generateFlippedOverlaps);
    }
}

int OverlapHifiWorkflow::Runner(const PacBio::CLI_v2::Results& options)
{
    OverlapHifiSettings settings{options};

    std::string targetSeqDBFile = settings.TargetDBPrefix + ".seqdb";
    std::string targetSeedDBFile = settings.TargetDBPrefix + ".seeddb";
    std::string querySeqDBFile = settings.QueryDBPrefix + ".seqdb";
    std::string querySeedDBFile = settings.QueryDBPrefix + ".seeddb";

    PBLOG_INFO << "Loading the input DBs.";
    TicToc ttInit;

    // Load the target DB caches.
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> targetSeqDBCache =
        PacBio::Pancake::LoadSeqDBIndexCache(targetSeqDBFile);
    PBLOG_INFO << "After loading target seq cache: " << ttInit.VerboseSecs(true);
    std::shared_ptr<PacBio::Pancake::SeedDBIndexCache> targetSeedDBCache =
        PacBio::Pancake::LoadSeedDBIndexCache(targetSeedDBFile);
    PBLOG_INFO << "After loading target seed cache: " << ttInit.VerboseSecs(true);
    PBLOG_INFO << "Target seed params: k = " << targetSeedDBCache->seedParams.KmerSize
               << ", w = " << targetSeedDBCache->seedParams.MinimizerWindow
               << ", s = " << targetSeedDBCache->seedParams.Spacing
               << ", hpc = " << targetSeedDBCache->seedParams.UseHPC
               << ", rc = " << targetSeedDBCache->seedParams.UseRC;
    if (targetSeedDBCache->seedParams.UseHPC != settings.UseHPC) {
        throw std::runtime_error(
            "The --use-hpc option was either used to compute the target SeedDB but not specified "
            "for overlapping, or vice versa.");
    }

    // Load the query DB caches.
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> querySeqDBCache =
        PacBio::Pancake::LoadSeqDBIndexCache(querySeqDBFile);
    PBLOG_INFO << "After loading query seq cache: " << ttInit.VerboseSecs(true);
    std::shared_ptr<PacBio::Pancake::SeedDBIndexCache> querySeedDBCache =
        PacBio::Pancake::LoadSeedDBIndexCache(querySeedDBFile);
    PBLOG_INFO << "After loading query seed cache: " << ttInit.VerboseSecs(true);
    PBLOG_INFO << "Query seed params: k = " << targetSeedDBCache->seedParams.KmerSize
               << ", w = " << targetSeedDBCache->seedParams.MinimizerWindow
               << ", s = " << targetSeedDBCache->seedParams.Spacing
               << ", hpc = " << targetSeedDBCache->seedParams.UseHPC
               << ", rc = " << targetSeedDBCache->seedParams.UseRC;
    if (querySeedDBCache->seedParams.UseHPC != settings.UseHPC) {
        throw std::runtime_error(
            "The --use-hpc option was either used to compute the query SeedDB but not specified "
            "for overlapping, or vice versa.");
    }

    // Create the target readers.
    PacBio::Pancake::SeqDBReaderCachedBlock targetSeqDBReader(targetSeqDBCache, settings.UseHPC);
    targetSeqDBReader.LoadBlocks({settings.TargetBlockId});
    PacBio::Pancake::SeedDBReaderRawBlock targetSeedDBReader(targetSeedDBCache);

    // Read the seeds for the target block.
    std::vector<PacBio::Pancake::SeedDB::SeedRaw> targetSeeds =
        targetSeedDBReader.GetBlock(settings.TargetBlockId);

    ttInit.Stop();
    PBLOG_INFO << "Loaded the target index and seqs in " << ttInit.GetSecs() << " sec / "
               << ttInit.GetCpuSecs() << " CPU sec";

    PBLOG_INFO << "Target seqs: " << targetSeqDBReader.records().size();
    PBLOG_INFO << "Target seeds: " << targetSeeds.size();

    // Build the seed index.
    TicToc ttIndex;
    PacBio::Pancake::SeedIndex index(targetSeedDBCache, std::move(targetSeeds));
    ttIndex.Stop();
    PBLOG_INFO << "Built the seed index in " << ttIndex.GetSecs() << " sec / "
               << ttIndex.GetCpuSecs() << " CPU sec";

    // Seed statistics, and computing the cutoff.
    TicToc ttSeedStats;
    int64_t freqMax = 0;
    int64_t freqCutoff = 0;
    double freqAvg = 0.0;
    double freqMedian = 0.0;
    index.ComputeFrequencyStats(settings.FreqPercentile, freqMax, freqAvg, freqMedian, freqCutoff);
    ttSeedStats.Stop();
    PBLOG_INFO << "Computed the seed frequency statistics in " << ttSeedStats.GetSecs() << " sec / "
               << ttSeedStats.GetCpuSecs() << " CPU sec";

    PBLOG_INFO << "Seed statistic: freqMax = " << freqMax << ", freqAvg = " << freqAvg
               << ", freqMedian = " << freqMedian << ", freqCutoff = " << freqCutoff;

    PBLOG_INFO << "Using traceback alignment: " << (settings.UseTraceback ? "on" : "off");

    PBLOG_INFO << "Beginning to map the sequences.";
    PBLOG_INFO << "Using " << settings.NumThreads << " threads.";
    std::vector<Mapper> mappers;
    for (size_t i = 0; i < settings.NumThreads; ++i) {
        mappers.emplace_back(Mapper(settings));
    }
    TicToc ttMap;

    auto writer = PacBio::Pancake::OverlapWriterFactory(settings.OutFormat, stdout,
                                                        settings.WriteIds, settings.WriteCigar);

    writer->WriteHeader(targetSeqDBReader);

    const int32_t endBlockId = (settings.QueryBlockEndId <= 0) ? querySeqDBCache->blockLines.size()
                                                               : settings.QueryBlockEndId;

    // Process all blocks.
    PacBio::Pancake::SeqDBReaderCachedBlock querySeqDBReader(querySeqDBCache, settings.UseHPC);
    PacBio::Pancake::SeedDBReaderCachedBlock querySeedDBReader(querySeedDBCache);
    for (int32_t queryBlockId = settings.QueryBlockStartId; queryBlockId < endBlockId;
         queryBlockId += settings.CombineBlocks) {

        std::vector<int32_t> blocksToLoad;
        for (int32_t blockId = queryBlockId;
             blockId < std::min(endBlockId, (queryBlockId + settings.CombineBlocks)); ++blockId) {
            blocksToLoad.emplace_back(blockId);
        }
        if (blocksToLoad.empty()) {
            throw std::runtime_error("There are zero blocks to load!");
        }
        std::string blocksToLoadStr;
        {
            std::ostringstream oss;
            oss << "{" << blocksToLoad[0];
            for (size_t i = 1; i < blocksToLoad.size(); ++i) {
                oss << ", " << blocksToLoad[i];
            }
            oss << "}";
            blocksToLoadStr = oss.str();
        }

        PBLOG_INFO << "Loading the query blocks: " << blocksToLoadStr << ".";
        // Create the query readers for the current block.
        TicToc ttQueryLoad;
        // PacBio::Pancake::SeqDBReaderCached querySeqDBReader(querySeqDBCache, queryBlockId);
        querySeqDBReader.LoadBlocks(blocksToLoad);
        PBLOG_INFO << "Loaded the query SeqDB cache block after " << ttQueryLoad.GetSecs(true)
                   << " sec / " << ttQueryLoad.GetCpuSecs(true) << " CPU sec";
        querySeedDBReader.LoadBlock(blocksToLoad);
        ttQueryLoad.Stop();
        PBLOG_INFO << "Loaded the query SeedDB cache block after " << ttQueryLoad.GetSecs()
                   << " sec / " << ttQueryLoad.GetCpuSecs() << " CPU sec";
        PBLOG_INFO << "Loaded all query blocks in " << ttQueryLoad.GetSecs() << " sec / "
                   << ttQueryLoad.GetCpuSecs() << " CPU sec";

        std::ostringstream oss;
        oss << "About to map query blocks: " << blocksToLoadStr
            << ": num_seqs = " << querySeqDBReader.records().size();
        PBLOG_INFO << oss.str();

        // Parallel processing.
        {
            TicToc ttQueryBlockMapping;
            const int32_t numRecords = static_cast<int32_t>(querySeqDBReader.records().size());
            // Storage for the results of the batch run.
            std::vector<PacBio::Pancake::MapperResult> results(numRecords);
            // No empty threads (e.g. reduce the requested number, if small record input)
            const int32_t actualThreadCount =
                std::min(static_cast<int32_t>(settings.NumThreads), numRecords);

            // Determine how many records should land in each thread, spread roughly evenly.
            const int32_t minimumRecordsPerThreads = (numRecords / actualThreadCount);
            const int32_t remainingRecords = (numRecords % actualThreadCount);
            std::vector<int32_t> recordsPerThread(actualThreadCount, minimumRecordsPerThreads);
            for (int i = 0; i < remainingRecords; ++i)
                ++recordsPerThread[i];

            // Run the mapping in parallel.
            PacBio::Parallel::FireAndForget faf(settings.NumThreads);
            int32_t submittedCount = 0;
            for (int32_t i = 0; i < actualThreadCount; ++i) {
                faf.ProduceWith(Worker, std::cref(targetSeqDBReader), std::cref(index),
                                std::cref(querySeqDBReader), std::cref(querySeedDBReader),
                                std::cref(settings), std::cref(mappers[i]), freqCutoff,
                                settings.WriteReverseOverlaps, submittedCount,
                                submittedCount + recordsPerThread[i], std::ref(results));
                submittedCount += recordsPerThread[i];
            }
            faf.Finalize();

            // Write the results.
            for (size_t i = 0; i < querySeqDBReader.records().size(); ++i) {
                const auto& result = results[i];
                const auto& querySeq = querySeqDBReader.records()[i];

                for (const auto& ovl : result.overlaps) {
                    writer->Write(ovl, targetSeqDBReader, querySeq, ovl->IsFlipped);
                }
            }

            ttQueryBlockMapping.Stop();

            PBLOG_INFO << "Mapped query block in " << ttQueryBlockMapping.GetSecs() << " sec / "
                       << ttQueryBlockMapping.GetCpuSecs() << " CPU sec.";
        }
    }
    ttMap.Stop();
    PBLOG_INFO << "Mapped all query blocks in " << ttMap.GetSecs() << " sec / "
               << ttMap.GetCpuSecs() << " CPU sec.";

    return EXIT_SUCCESS;
}

}  // namespace Pancake
}  // namespace PacBio
