// Authors: Ivan Sovic

#include "seeddb/SeedDBWorkflow.h"
#include <seeddb/Minimizers.h>
#include <seeddb/Seed.h>
#include <seeddb/SeedDBWriter.h>
#include <seqdb/FastaSequenceId.h>
#include <seqdb/SeqDBIndexCache.h>
#include <seqdb/SeqDBReader.h>
#include "seeddb/SeedDBSettings.h"

#include <pbcopper/parallel/FireAndForget.h>
#include <pbcopper/parallel/WorkQueue.h>

#include <iostream>

namespace PacBio {
namespace Pancake {
namespace SeedDB {

void Worker(const std::vector<PacBio::Pancake::FastaSequenceId>& records,
            const SeedDBSettings& settings, int32_t start, int32_t end, int32_t startAbs,
            std::vector<std::vector<__int128>>& seeds)
{
    const auto& sp = settings.SeedParameters;

    for (int32_t i = start; i < end; ++i) {
        const auto& record = records[i];
        const uint8_t* seq = reinterpret_cast<const uint8_t*>(record.Bases().data());
        int32_t seqLen = record.Bases().size();
        int rv =
            GenerateMinimizers(seeds[i], seq, seqLen, 0, record.Id(), sp.KmerSize,
                               sp.MinimizerWindow, sp.Spacing, sp.UseRC, sp.UseHPC, sp.MaxHPCLen);
        if (rv)
            throw std::runtime_error("Generating minimizers failed, startAbs = " +
                                     std::to_string(startAbs));
    }
}

int SeedDBWorkflow::Runner(const PacBio::CLI_v2::Results& options)
{
    SeedDBSettings settings{options};

    // Load the DB.
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBCache =
        PacBio::Pancake::LoadSeqDBIndexCache(settings.InputFile);

    // Create a reader.
    PacBio::Pancake::SeqDBReader reader(seqDBCache);

    int32_t numBlocks = seqDBCache->blockLines.size();
    int32_t absOffset = 0;

    auto writer = PacBio::Pancake::CreateSeedDBWriter(settings.OutputPrefix, settings.SplitBlocks,
                                                      settings.SeedParameters);

    for (int32_t blockId = 0; blockId < numBlocks; ++blockId) {
        // Load a block of records.
        std::vector<PacBio::Pancake::FastaSequenceId> records;
        bool rv = reader.GetBlock(records, blockId);
        if (rv == false)
            throw std::runtime_error("Something went wrong when fetching block " +
                                     std::to_string(blockId));
        int32_t numRecords = records.size();

        // Generate seeds in parallel.
        std::vector<std::vector<__int128>> results(numRecords);
        PacBio::Parallel::FireAndForget faf(settings.NumThreads);
        for (int32_t i = 0; i < numRecords; ++i) {
            faf.ProduceWith(Worker, std::cref(records), std::cref(settings), i, i + 1,
                            i + absOffset, std::ref(results));
        }
        faf.Finalize();

        writer->WriteSeeds(records, results);
        writer->MarkBlockEnd();

        // Increase the abs offset counter.
        absOffset += numRecords;
    }

    return EXIT_SUCCESS;
}

}  // namespace SeedDB
}  // namespace Pancake
}  // namespace PacBio
