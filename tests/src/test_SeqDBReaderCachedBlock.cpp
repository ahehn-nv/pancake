// Authors: Ivan Sovic

#include <PancakeTestData.h>
#include <gtest/gtest.h>
#include <pacbio/seqdb/SeqDBReader.h>
#include <pacbio/seqdb/SeqDBReaderCachedBlock.h>
#include <iostream>

TEST(SeqDBReaderCachedBlock, BatchCompareWithSeqDBReader_UncompressedInput)
{
    /*
     * Batch test loading several SeqDBs in a row.
     * This test compares the same DB loaded with the Unit-Under-Test and with
     * a previously tested SeqDBReader. The loaded sequences should be the same.
     * The only difference is the format in which they return the sequences,
     * so they need to be collected in a comparable format first.
    */

    const std::vector<std::string> inDBs = {
        PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/test-3.seqdb",
        PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/test-7-uncompressed-2blocks.seqdb",
    };

    for (const auto& inSeqDB : inDBs) {
        std::cerr << "Testing DB: '" << inSeqDB << "'.\n";

        // Load the SeedDB.
        std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBCache =
            PacBio::Pancake::LoadSeqDBIndexCache(inSeqDB);

        const int32_t numBlocks = seqDBCache->blockLines.size();

        for (int32_t blockId = 0; blockId < numBlocks; ++blockId) {
            std::cerr << "  Evaluating block " << blockId << ".\n";

            // "Reference" (or "truth") reader. This was tested earlier, thoroughly.
            // Collect the expected results for this block using a trusty reader.
            PacBio::Pancake::SeqDBReader readerTruth(seqDBCache);
            std::vector<PacBio::Pancake::FastaSequenceId> expected;
            readerTruth.GetBlock(expected, blockId);

            // Create a SeedDB reader under test.
            // Convert all FastaSequenceCached to FastaSequenceId for easier comparison.
            PacBio::Pancake::SeqDBReaderCachedBlock readerTest(seqDBCache, blockId);
            std::vector<PacBio::Pancake::FastaSequenceId> results;
            for (const auto& record : readerTest.records()) {
                results.emplace_back(PacBio::Pancake::FastaSequenceId(
                    record.Name(), std::string(record.Bases(), record.Size()), record.Id()));
            }

            // Evaluate the current block.
            EXPECT_EQ(expected, results);
        }
    }
}

TEST(SeqDBReaderCachedBlock, BatchCompareWithSeqDBReader_CompressedInput)
{
    /*
     * Same as before, but the input DBs are compressed.
    */

    const std::vector<std::string> inDBs = {
        PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/test-1.seqdb",
        PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/test-6.seqdb",
        PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/test-8-compressed-2blocks.seqdb",
    };

    for (const auto& inSeqDB : inDBs) {
        std::cerr << "Testing DB: '" << inSeqDB << "'.\n";

        // Load the SeedDB.
        std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBCache =
            PacBio::Pancake::LoadSeqDBIndexCache(inSeqDB);

        const int32_t numBlocks = seqDBCache->blockLines.size();

        for (int32_t blockId = 0; blockId < numBlocks; ++blockId) {
            std::cerr << "  Evaluating block " << blockId << ".\n";

            // "Reference" (or "truth") reader. This was tested earlier, thoroughly.
            // Collect the expected results for this block using a trusty reader.
            PacBio::Pancake::SeqDBReader readerTruth(seqDBCache);
            std::vector<PacBio::Pancake::FastaSequenceId> expected;
            readerTruth.GetBlock(expected, blockId);

            // Create a SeedDB reader under test.
            // Convert all FastaSequenceCached to FastaSequenceId for easier comparison.
            PacBio::Pancake::SeqDBReaderCachedBlock readerTest(seqDBCache, blockId);
            std::vector<PacBio::Pancake::FastaSequenceId> results;
            for (const auto& record : readerTest.records()) {
                results.emplace_back(PacBio::Pancake::FastaSequenceId(
                    record.Name(), std::string(record.Bases(), record.Size()), record.Id()));
            }

            // Evaluate the current block.
            EXPECT_EQ(expected, results);
        }
    }
}
