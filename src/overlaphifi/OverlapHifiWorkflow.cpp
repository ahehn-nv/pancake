// Authors: Ivan Sovic

#include "overlaphifi/OverlapHifiWorkflow.h"
#include <pacbio/overlaphifi/Mapper.h>
#include <pacbio/overlaphifi/OverlapHifiSettings.h>
#include <pacbio/overlaphifi/OverlapWriter.h>
#include <pacbio/overlaphifi/SeedIndex.h>
#include <pacbio/seeddb/Seed.h>
#include <pacbio/seeddb/SeedDBIndexCache.h>
#include <pacbio/seeddb/SeedDBReaderCached.h>
#include <pacbio/seeddb/SeedDBReaderRawBlock.h>
#include <pacbio/seqdb/SeqDBIndexCache.h>
#include <pacbio/seqdb/SeqDBReaderCached.h>
#include <pacbio/util/TicToc.h>
#include <pbcopper/logging/LogLevel.h>
#include <pbcopper/logging/Logging.h>
#include <pbcopper/parallel/FireAndForget.h>
#include <pbcopper/parallel/WorkQueue.h>

namespace PacBio {
namespace Pancake {

void Worker(const PacBio::Pancake::SeqDBReaderCached& targetSeqDBReader,
            const PacBio::Pancake::SeedIndex& index,
            const PacBio::Pancake::SeqDBReaderCached& querySeqDBReader,
            const PacBio::Pancake::SeedDBReaderCached& querySeedDBReader,
            const OverlapHifiSettings& /*settings*/, const PacBio::Pancake::Mapper& mapper,
            int64_t freqCutoff, int32_t start, int32_t end,
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
        results[i] = mapper.Map(targetSeqDBReader, index, querySeq, querySeeds, freqCutoff);
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
    std::shared_ptr<PacBio::Pancake::SeedDBIndexCache> targetSeedDBCache =
        PacBio::Pancake::LoadSeedDBIndexCache(targetSeedDBFile);
    // Load the query DB caches.
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> querySeqDBCache =
        PacBio::Pancake::LoadSeqDBIndexCache(querySeqDBFile);
    std::shared_ptr<PacBio::Pancake::SeedDBIndexCache> querySeedDBCache =
        PacBio::Pancake::LoadSeedDBIndexCache(querySeedDBFile);
    // Create the target readers.
    PacBio::Pancake::SeqDBReaderCached targetSeqDBReader(targetSeqDBCache, settings.TargetBlockId);
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

    PBLOG_INFO << "Beginning to map the sequences.";
    PBLOG_INFO << "Using " << settings.NumThreads << " threads.";
    Mapper mapper(settings);
    TicToc ttMap;

    PacBio::Pancake::OverlapWriter writer(stdout, settings.WriteReverseOverlaps, settings.WriteIds);

    const int32_t endBlockId = (settings.QueryBlockEndId <= 0) ? querySeqDBCache->blockLines.size()
                                                               : settings.QueryBlockEndId;

    // Process all blocks.
    for (int32_t queryBlockId = settings.QueryBlockStartId; queryBlockId < endBlockId;
         ++queryBlockId) {
        PBLOG_INFO << "Loading the query block " << queryBlockId << ".";
        // Create the query readers for the current block.
        TicToc ttQueryLoad;
        PacBio::Pancake::SeqDBReaderCached querySeqDBReader(querySeqDBCache, queryBlockId);
        PacBio::Pancake::SeedDBReaderCached querySeedDBReader(querySeedDBCache, queryBlockId);
        ttQueryLoad.Stop();
        PBLOG_INFO << "Loaded the query block in " << ttQueryLoad.GetSecs() << " sec / "
                   << ttQueryLoad.GetCpuSecs() << " CPU sec";

        PBLOG_INFO << "About to map query block " << queryBlockId
                   << ": num_seqs = " << querySeqDBReader.records().size();
        // Parallel processing.
        {
            TicToc ttQueryBlockMapping;
            const int32_t numRecords = static_cast<int32_t>(querySeqDBReader.records().size());
            std::vector<PacBio::Pancake::MapperResult> results(numRecords);

            // Run the mapping in parallel.
            PacBio::Parallel::FireAndForget faf(settings.NumThreads);
            for (int32_t i = 0; i < numRecords; ++i) {
                faf.ProduceWith(Worker, std::cref(targetSeqDBReader), std::cref(index),
                                std::cref(querySeqDBReader), std::cref(querySeedDBReader),
                                std::cref(settings), std::cref(mapper), freqCutoff, i, i + 1,
                                std::ref(results));
            }
            faf.Finalize();

            // Write the results.
            for (size_t i = 0; i < querySeqDBReader.records().size(); ++i) {
                const auto& result = results[i];
                const auto& querySeq = querySeqDBReader.records()[i];
                writer.Write(result.overlaps, targetSeqDBReader, querySeq);
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
