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
        // Load the SeedDB.
        std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBCache =
            PacBio::Pancake::LoadSeqDBIndexCache(inSeqDB);

        const int32_t numBlocks = seqDBCache->blockLines.size();

        for (int32_t blockId = 0; blockId < numBlocks; ++blockId) {
            SCOPED_TRACE(inSeqDB + ", blockId = " + std::to_string(blockId));

            // "Reference" (or "truth") reader. This was tested earlier, thoroughly.
            // Collect the expected results for this block using a trusty reader.
            PacBio::Pancake::SeqDBReader readerTruth(seqDBCache);
            std::vector<PacBio::Pancake::FastaSequenceId> expected;
            readerTruth.GetBlock(expected, blockId);

            // Create a SeedDB reader under test.
            // Convert all FastaSequenceCached to FastaSequenceId for easier comparison.
            PacBio::Pancake::SeqDBReaderCachedBlock readerTest(seqDBCache, {blockId});
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
        // Load the SeedDB.
        std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBCache =
            PacBio::Pancake::LoadSeqDBIndexCache(inSeqDB);

        const int32_t numBlocks = seqDBCache->blockLines.size();

        for (int32_t blockId = 0; blockId < numBlocks; ++blockId) {
            SCOPED_TRACE(inSeqDB + ", blockId = " + std::to_string(blockId));

            // "Reference" (or "truth") reader. This was tested earlier, thoroughly.
            // Collect the expected results for this block using a trusty reader.
            PacBio::Pancake::SeqDBReader readerTruth(seqDBCache);
            std::vector<PacBio::Pancake::FastaSequenceId> expected;
            readerTruth.GetBlock(expected, blockId);

            // Create a SeedDB reader under test.
            // Convert all FastaSequenceCached to FastaSequenceId for easier comparison.
            PacBio::Pancake::SeqDBReaderCachedBlock readerTest(seqDBCache, {blockId});
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

TEST(SeqDBReaderCachedBlock, MultipleInputBlocks_Uncompressed)
{
    /*
     * Same as before, but the input DBs are compressed.
    */

    const std::string inSeqDB = PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/test-3.seqdb";
    const std::vector<int32_t> inBlocks = {1, 2, 3};

    // Load the SeedDB.
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBCache =
        PacBio::Pancake::LoadSeqDBIndexCache(inSeqDB);

    // Collect all expected sequences for the specified input blocks
    // using an orthogonal reader. These are treated as the truth sequences.
    std::vector<PacBio::Pancake::FastaSequenceId> expected;
    for (const auto& blockId : inBlocks) {
        // "Reference" (or "truth") reader. This was tested earlier, thoroughly.
        // Collect the expected results for this block using a trusty reader.
        std::vector<PacBio::Pancake::FastaSequenceId> currExpected;
        PacBio::Pancake::SeqDBReader readerTruth(seqDBCache);
        readerTruth.GetBlock(currExpected, blockId);

        for (const auto& val : currExpected) {
            expected.emplace_back(val);
        }
    }

    // Create a unit under test.
    // Read the sequences for the specified blocks, and convert all
    // FastaSequenceCached to FastaSequenceId for easier comparison.
    PacBio::Pancake::SeqDBReaderCachedBlock readerTest(seqDBCache, inBlocks);
    std::vector<PacBio::Pancake::FastaSequenceId> results;
    for (const auto& record : readerTest.records()) {
        results.emplace_back(PacBio::Pancake::FastaSequenceId(
            record.Name(), std::string(record.Bases(), record.Size()), record.Id()));
    }

    // Evaluate the current block.
    EXPECT_EQ(expected, results);
}

TEST(SeqDBReaderCachedBlock, MultipleInputBlocks_Compressed)
{
    /*
     * Same as before, but the input DBs are compressed.
    */

    const std::string inSeqDB = PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/test-1.seqdb";
    const std::vector<int32_t> inBlocks = {1, 2, 3};

    // Load the SeedDB.
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBCache =
        PacBio::Pancake::LoadSeqDBIndexCache(inSeqDB);

    // Collect all expected sequences with an orthogonal reader.
    std::vector<PacBio::Pancake::FastaSequenceId> expected;
    PacBio::Pancake::SeqDBReader readerTruth(seqDBCache);
    for (const auto& blockId : inBlocks) {
        // "Reference" (or "truth") reader. This was tested earlier, thoroughly.
        // Collect the expected results for this block using a trusty reader.
        std::vector<PacBio::Pancake::FastaSequenceId> currExpected;
        readerTruth.GetBlock(currExpected, blockId);

        for (const auto& val : currExpected) {
            expected.emplace_back(val);
        }
    }

    // Create a unit under test.
    // Convert all FastaSequenceCached to FastaSequenceId for easier comparison.
    PacBio::Pancake::SeqDBReaderCachedBlock readerTest(seqDBCache, inBlocks);
    std::vector<PacBio::Pancake::FastaSequenceId> results;
    for (const auto& record : readerTest.records()) {
        results.emplace_back(PacBio::Pancake::FastaSequenceId(
            record.Name(), std::string(record.Bases(), record.Size()), record.Id()));
    }

    // Evaluate the current block.
    EXPECT_EQ(expected, results);
}
