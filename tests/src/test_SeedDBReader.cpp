// Authors: Ivan Sovic

#include <PancakeTestData.h>
#include <gtest/gtest.h>
#include <pacbio/seeddb/Minimizers.h>
#include <pacbio/seeddb/SeedDBReader.h>
#include <pacbio/seeddb/SequenceSeeds.h>
#include <pacbio/seqdb/SeqDBReader.h>
#include <seqdb/FastaSequenceId.h>
#include <sstream>

// std::vector<PacBio::BAM::FastaSequence> HelperLoadFasta(const std::string& inFasta)
// {
//     std::vector<PacBio::BAM::FastaSequence> ret;
//     PacBio::BAM::FastaReader inReader{inFasta};
//     PacBio::BAM::FastaSequence record;
//     while (inReader.GetNext(record))
//         ret.emplace_back(record);
//     return ret;
// }

// std::string HelperLoadFastaAsString(const std::string& inFasta)
// {
//     std::ostringstream oss;
//     auto records = HelperLoadFasta(inFasta);
//     for (const auto& record : records)
//         oss << ">" << record.Name() << "\n" << record.Bases() << "\n";
//     return oss.str();
// }

std::string HelperLoadFileAsString(const std::string& in_path)
{
    std::ifstream ifs(in_path);
    std::string line;
    std::ostringstream oss;
    while (std::getline(ifs, line)) {
        oss << line << "\n";
    }
    return oss.str();
}

std::vector<PacBio::Pancake::FastaSequenceId> HelperLoadAllFromSeqDB(const std::string& inSeqDB)
{
    // Load the DB.
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBCache =
        PacBio::Pancake::LoadSeqDBIndexCache(inSeqDB);

    // Create a reader.
    PacBio::Pancake::SeqDBReader reader(seqDBCache);

    // Load all sequences.
    std::vector<PacBio::Pancake::FastaSequenceId> records;
    reader.GetNextBatch(records, -1);
    return records;
}

PacBio::Pancake::SequenceSeeds HelperComputeSequenceSeeds(
    const PacBio::Pancake::FastaSequenceId& record, int32_t k, int32_t w, bool useHPC,
    int32_t maxHPCLen, bool useRC)
{
    const uint8_t* seqData = reinterpret_cast<const uint8_t*>(record.Bases().data());
    int32_t seqLen = record.Bases().size();
    std::vector<__int128> seeds;
    PacBio::Pancake::SeedDB::GenerateMinimizers(seeds, seqData, seqLen, 0, record.Id(), k, w, useRC,
                                                useHPC, maxHPCLen);
    PacBio::Pancake::SequenceSeeds seqSeeds(record.Name(), std::move(seeds), record.Id());
    return seqSeeds;
}

std::vector<PacBio::Pancake::SequenceSeeds> HelperComputeSequenceSeeds(
    const std::vector<PacBio::Pancake::FastaSequenceId>& records, int32_t k, int32_t w, bool useHPC,
    int32_t maxHPCLen, bool useRC)
{
    std::vector<PacBio::Pancake::SequenceSeeds> ret;
    for (const auto& record : records) {
        ret.emplace_back(HelperComputeSequenceSeeds(record, k, w, useHPC, maxHPCLen, useRC));
    }
    return ret;
}

TEST(SeedDBIndexCache, RoundTrip1)
{
    /*
     * This test parses the SeedDB index and dumps it into a string.
     * The string shuld match the exact text of the SeedDB index file.
    */

    // Input values;
    const std::string inSeedDB =
        PacBio::PancakeTestsConfig::Data_Dir + "/seeddb-writer/test-1a.seeddb";

    // Expected output.
    const auto expected = HelperLoadFileAsString(inSeedDB);

    // Load the SeedDB.
    std::shared_ptr<PacBio::Pancake::SeedDBIndexCache> seedDBCache =
        PacBio::Pancake::LoadSeedDBIndexCache(inSeedDB);

    std::ostringstream result;
    result << *seedDBCache;

    EXPECT_EQ(expected, result.str());
}

TEST(SeedDBIndexCache, RoundTrip2)
{
    /*
     * Same as RoundTrip1 but on a different input DB.
     * This test parses the SeedDB index and dumps it into a string.
     * The string shuld match the exact text of the SeedDB index file.
    */

    // Input values;
    const std::string inSeedDB =
        PacBio::PancakeTestsConfig::Data_Dir + "/seeddb-writer/test-1b.seeddb";

    // Expected output.
    const auto expected = HelperLoadFileAsString(inSeedDB);

    // Load the SeedDB.
    std::shared_ptr<PacBio::Pancake::SeedDBIndexCache> seedDBCache =
        PacBio::Pancake::LoadSeedDBIndexCache(inSeedDB);

    std::ostringstream result;
    result << *seedDBCache;

    EXPECT_EQ(expected, result.str());
}

TEST(SeedDBIndexCache, NonexistentDB)
{
    /*
     * Same as RoundTrip1 but on a different input DB.
     * This test parses the SeedDB index and dumps it into a string.
     * The string shuld match the exact text of the SeedDB index file.
    */

    // Input values.
    const std::string inSeedDB = PacBio::PancakeTestsConfig::Data_Dir + "/nonexistent.seeddb";

    EXPECT_THROW(
        {
            // Load the DB.
            std::shared_ptr<PacBio::Pancake::SeedDBIndexCache> seedDBCache =
                PacBio::Pancake::LoadSeedDBIndexCache(inSeedDB);
        },
        std::runtime_error);
}

TEST(SeedDBReader, GetNext1)
{
    /*
     * Iterates through all records in the SeedDB and loads the SequenceSeeds
     * records one by one.
     *
     * Expected results - This test manually fetches sequences from a
     * SeqDB and computes the seeds for those sequences. The computed seeds
     * should be the same as the loaded ones.
     *
     * The GenerateMinimizers function is tested elsewhere.
     * The SeqDBReader::GetNext is also tested elsewhere.
    */

    // Input values.
    const std::string inSeqDB = PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/test-1.seqdb";
    const std::string inSeedDB =
        PacBio::PancakeTestsConfig::Data_Dir + "/seeddb-writer/test-1a.seeddb";

    // Load the SeedDB.
    std::shared_ptr<PacBio::Pancake::SeedDBIndexCache> seedDBCache =
        PacBio::Pancake::LoadSeedDBIndexCache(inSeedDB);

    // Expected results - compute the minimizers from the sequences here.
    std::vector<PacBio::Pancake::FastaSequenceId> seqDBRecords = HelperLoadAllFromSeqDB(inSeqDB);
    std::vector<PacBio::Pancake::SequenceSeeds> expected = HelperComputeSequenceSeeds(
        seqDBRecords, seedDBCache->seedParams.KmerSize, seedDBCache->seedParams.MinimizerWindow,
        seedDBCache->seedParams.UseHPC, seedDBCache->seedParams.MaxHPCLen,
        seedDBCache->seedParams.UseRC);

    // Create a SeedDB reader.
    PacBio::Pancake::SeedDBReader reader(seedDBCache);

    // Collect all seeds from the SeedDB.
    std::vector<PacBio::Pancake::SequenceSeeds> results;
    PacBio::Pancake::SequenceSeeds record;
    while (reader.GetNext(record)) {
        results.emplace_back(record);
    }

    // Evaluate.
    EXPECT_EQ(expected, results);
}

TEST(SeedDBReader, GetNext2)
{
    /*
     * Same as before, but tests loading seeds from multiple .seeds files.
     *
     * Iterates through all records in the SeedDB and loads the SequenceSeeds
     * records one by one.
     *
     * Expected results - This test manually fetches sequences from a
     * SeqDB and computes the seeds for those sequences. The computed seeds
     * should be the same as the loaded ones.
     *
     * The GenerateMinimizers function is tested elsewhere.
     * The SeqDBReader::GetNext is also tested elsewhere.
    */

    // Input values.
    const std::string inSeqDB = PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/test-1.seqdb";
    const std::string inSeedDB =
        PacBio::PancakeTestsConfig::Data_Dir + "/seeddb-writer/test-1b.seeddb";

    // Load the SeedDB.
    std::shared_ptr<PacBio::Pancake::SeedDBIndexCache> seedDBCache =
        PacBio::Pancake::LoadSeedDBIndexCache(inSeedDB);

    // Expected results - compute the minimizers from the sequences here.
    std::vector<PacBio::Pancake::FastaSequenceId> seqDBRecords = HelperLoadAllFromSeqDB(inSeqDB);
    std::vector<PacBio::Pancake::SequenceSeeds> expected = HelperComputeSequenceSeeds(
        seqDBRecords, seedDBCache->seedParams.KmerSize, seedDBCache->seedParams.MinimizerWindow,
        seedDBCache->seedParams.UseHPC, seedDBCache->seedParams.MaxHPCLen,
        seedDBCache->seedParams.UseRC);

    // Create a SeedDB reader.
    PacBio::Pancake::SeedDBReader reader(seedDBCache);

    // Collect all seeds from the SeedDB.
    std::vector<PacBio::Pancake::SequenceSeeds> results;
    PacBio::Pancake::SequenceSeeds record;
    while (reader.GetNext(record)) {
        results.emplace_back(record);
    }

    // Evaluate.
    EXPECT_EQ(expected, results);
}

TEST(SeedDBReader, GetSeedsByID1)
{
    /*
     * Iterate through expected sequences in reverse, and see if the loaded
     * is the same as expected.
     *
     * Expected results - This test manually fetches sequences from a
     * SeqDB and computes the seeds for those sequences. The computed seeds
     * should be the same as the loaded ones.
     *
     * The GenerateMinimizers function is tested elsewhere.
     * The SeqDBReader::GetNext is also tested elsewhere.
    */

    // Input values.
    const std::string inSeqDB = PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/test-1.seqdb";
    const std::string inSeedDB =
        PacBio::PancakeTestsConfig::Data_Dir + "/seeddb-writer/test-1a.seeddb";

    // Load the SeedDB.
    std::shared_ptr<PacBio::Pancake::SeedDBIndexCache> seedDBCache =
        PacBio::Pancake::LoadSeedDBIndexCache(inSeedDB);

    // Expected results - compute the minimizers from the sequences here.
    std::vector<PacBio::Pancake::FastaSequenceId> seqDBRecords = HelperLoadAllFromSeqDB(inSeqDB);
    std::vector<PacBio::Pancake::SequenceSeeds> expected = HelperComputeSequenceSeeds(
        seqDBRecords, seedDBCache->seedParams.KmerSize, seedDBCache->seedParams.MinimizerWindow,
        seedDBCache->seedParams.UseHPC, seedDBCache->seedParams.MaxHPCLen,
        seedDBCache->seedParams.UseRC);

    // Create a SeedDB reader.
    PacBio::Pancake::SeedDBReader reader(seedDBCache);

    // Reverse access all SequenceSeeds and fetch the ID from the SeedDB.
    // Check if they're the same.
    for (int32_t i = static_cast<int32_t>(expected.size()) - 1; i >= 0; --i) {
        const auto& expRecord = expected[i];

        PacBio::Pancake::SequenceSeeds record;
        reader.GetSeedsForSequence(record, expRecord.Id());

        EXPECT_EQ(expRecord, record);
    }
}

TEST(SeedDBReader, GetSeedsByID2)
{
    /*
     * Same as before, but tests loading seeds from multiple .seeds files.
     *
     * Iterate through expected sequences in reverse, and see if the loaded
     * is the same as expected.
     *
     * Expected results - This test manually fetches sequences from a
     * SeqDB and computes the seeds for those sequences. The computed seeds
     * should be the same as the loaded ones.
     *
     * The GenerateMinimizers function is tested elsewhere.
     * The SeqDBReader::GetNext is also tested elsewhere.
    */

    // Input values.
    const std::string inSeqDB = PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/test-1.seqdb";
    const std::string inSeedDB =
        PacBio::PancakeTestsConfig::Data_Dir + "/seeddb-writer/test-1b.seeddb";

    // Load the SeedDB.
    std::shared_ptr<PacBio::Pancake::SeedDBIndexCache> seedDBCache =
        PacBio::Pancake::LoadSeedDBIndexCache(inSeedDB);

    // Expected results - compute the minimizers from the sequences here.
    std::vector<PacBio::Pancake::FastaSequenceId> seqDBRecords = HelperLoadAllFromSeqDB(inSeqDB);
    std::vector<PacBio::Pancake::SequenceSeeds> expected = HelperComputeSequenceSeeds(
        seqDBRecords, seedDBCache->seedParams.KmerSize, seedDBCache->seedParams.MinimizerWindow,
        seedDBCache->seedParams.UseHPC, seedDBCache->seedParams.MaxHPCLen,
        seedDBCache->seedParams.UseRC);

    // Create a SeedDB reader.
    PacBio::Pancake::SeedDBReader reader(seedDBCache);

    // Reverse access all SequenceSeeds and fetch the ID from the SeedDB.
    // Check if they're the same.
    for (int32_t i = static_cast<int32_t>(expected.size()) - 1; i >= 0; --i) {
        const auto& expRecord = expected[i];

        PacBio::Pancake::SequenceSeeds record;
        reader.GetSeedsForSequence(record, expRecord.Id());

        EXPECT_EQ(expRecord, record);
    }
}

TEST(SeedDBReader, GetSeedsByName1)
{
    /*
     * Iterate through expected sequences in reverse, and see if the loaded
     * is the same as expected.
     *
     * Expected results - This test manually fetches sequences from a
     * SeqDB and computes the seeds for those sequences. The computed seeds
     * should be the same as the loaded ones.
     *
     * The GenerateMinimizers function is tested elsewhere.
     * The SeqDBReader::GetNext is also tested elsewhere.
    */

    // Input values.
    const std::string inSeqDB = PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/test-1.seqdb";
    const std::string inSeedDB =
        PacBio::PancakeTestsConfig::Data_Dir + "/seeddb-writer/test-1a.seeddb";

    // Load the SeedDB.
    std::shared_ptr<PacBio::Pancake::SeedDBIndexCache> seedDBCache =
        PacBio::Pancake::LoadSeedDBIndexCache(inSeedDB);

    // Expected results - compute the minimizers from the sequences here.
    std::vector<PacBio::Pancake::FastaSequenceId> seqDBRecords = HelperLoadAllFromSeqDB(inSeqDB);
    std::vector<PacBio::Pancake::SequenceSeeds> expected = HelperComputeSequenceSeeds(
        seqDBRecords, seedDBCache->seedParams.KmerSize, seedDBCache->seedParams.MinimizerWindow,
        seedDBCache->seedParams.UseHPC, seedDBCache->seedParams.MaxHPCLen,
        seedDBCache->seedParams.UseRC);

    // Create a SeedDB reader.
    PacBio::Pancake::SeedDBReader reader(seedDBCache);

    // Reverse access all SequenceSeeds and fetch the ID from the SeedDB.
    // Check if they're the same.
    for (int32_t i = static_cast<int32_t>(expected.size()) - 1; i >= 0; --i) {
        const auto& expRecord = expected[i];

        PacBio::Pancake::SequenceSeeds record;
        reader.GetSeedsForSequence(record, expRecord.Name());

        EXPECT_EQ(expRecord, record);
    }
}

TEST(SeedDBReader, GetSeedsByName2)
{
    /*
     * Same as before, but tests loading seeds from multiple .seeds files.
     *
     * Iterate through expected sequences in reverse, and see if the loaded
     * is the same as expected.
     *
     * Expected results - This test manually fetches sequences from a
     * SeqDB and computes the seeds for those sequences. The computed seeds
     * should be the same as the loaded ones.
     *
     * The GenerateMinimizers function is tested elsewhere.
     * The SeqDBReader::GetNext is also tested elsewhere.
    */

    // Input values.
    const std::string inSeqDB = PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/test-1.seqdb";
    const std::string inSeedDB =
        PacBio::PancakeTestsConfig::Data_Dir + "/seeddb-writer/test-1b.seeddb";

    // Load the SeedDB.
    std::shared_ptr<PacBio::Pancake::SeedDBIndexCache> seedDBCache =
        PacBio::Pancake::LoadSeedDBIndexCache(inSeedDB);

    // Expected results - compute the minimizers from the sequences here.
    std::vector<PacBio::Pancake::FastaSequenceId> seqDBRecords = HelperLoadAllFromSeqDB(inSeqDB);
    std::vector<PacBio::Pancake::SequenceSeeds> expected = HelperComputeSequenceSeeds(
        seqDBRecords, seedDBCache->seedParams.KmerSize, seedDBCache->seedParams.MinimizerWindow,
        seedDBCache->seedParams.UseHPC, seedDBCache->seedParams.MaxHPCLen,
        seedDBCache->seedParams.UseRC);

    // Create a SeedDB reader.
    PacBio::Pancake::SeedDBReader reader(seedDBCache);

    // Reverse access all SequenceSeeds and fetch the ID from the SeedDB.
    // Check if they're the same.
    for (int32_t i = static_cast<int32_t>(expected.size()) - 1; i >= 0; --i) {
        const auto& expRecord = expected[i];

        PacBio::Pancake::SequenceSeeds record;
        reader.GetSeedsForSequence(record, expRecord.Name());

        EXPECT_EQ(expRecord, record);
    }
}

TEST(SeedDBReader, GetSeedsNotExisting)
{
    /*
     * Try loading a record for a sequence name which does not exist in the SeedDB.
     * This should throw.
     *
     * The GenerateMinimizers function is tested elsewhere.
     * The SeqDBReader::GetNext is also tested elsewhere.
    */

    // Input values.
    const std::string inSeedDB =
        PacBio::PancakeTestsConfig::Data_Dir + "/seeddb-writer/test-1a.seeddb";

    // Load the SeedDB.
    std::shared_ptr<PacBio::Pancake::SeedDBIndexCache> seedDBCache =
        PacBio::Pancake::LoadSeedDBIndexCache(inSeedDB);

    // Create a SeedDB reader.
    PacBio::Pancake::SeedDBReader reader(seedDBCache);

    // Fetching non-existent sequence names throws.
    EXPECT_THROW(
        {
            PacBio::Pancake::SequenceSeeds record;
            reader.GetSeedsForSequence(record, "nonexistent-name");
        },
        std::runtime_error);

    // Fetching non-existent sequence ID throws.
    EXPECT_THROW(
        {
            PacBio::Pancake::SequenceSeeds record;
            reader.GetSeedsForSequence(record, 123);
        },
        std::runtime_error);
    EXPECT_THROW(
        {
            PacBio::Pancake::SequenceSeeds record;
            reader.GetSeedsForSequence(record, -1);
        },
        std::runtime_error);
}

TEST(SeedDBReader, GetNextBatch1)
{
    /*
     * Loads all records in one batch (batchSize == -1).
     *
     * Expected results - This test manually fetches sequences from a
     * SeqDB and computes the seeds for those sequences. The computed seeds
     * should be the same as the loaded ones.
     *
     * The GenerateMinimizers function is tested elsewhere.
     * The SeqDBReader::GetNext is also tested elsewhere.
    */

    // Input values.
    const std::string inSeqDB = PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/test-1.seqdb";
    const std::string inSeedDB =
        PacBio::PancakeTestsConfig::Data_Dir + "/seeddb-writer/test-1a.seeddb";

    // Load the SeedDB.
    std::shared_ptr<PacBio::Pancake::SeedDBIndexCache> seedDBCache =
        PacBio::Pancake::LoadSeedDBIndexCache(inSeedDB);

    // Expected results - compute the minimizers from the sequences here.
    std::vector<PacBio::Pancake::FastaSequenceId> seqDBRecords = HelperLoadAllFromSeqDB(inSeqDB);
    std::vector<PacBio::Pancake::SequenceSeeds> expected = HelperComputeSequenceSeeds(
        seqDBRecords, seedDBCache->seedParams.KmerSize, seedDBCache->seedParams.MinimizerWindow,
        seedDBCache->seedParams.UseHPC, seedDBCache->seedParams.MaxHPCLen,
        seedDBCache->seedParams.UseRC);

    // Create a SeedDB reader.
    PacBio::Pancake::SeedDBReader reader(seedDBCache);

    // Collect all seeds from the SeedDB.
    std::vector<PacBio::Pancake::SequenceSeeds> results;
    reader.GetNextBatch(results, -1);

    // Evaluate.
    EXPECT_EQ(expected, results);
}

TEST(SeedDBReader, GetNextBatch2)
{
    /*
     * Loads all batches of 10kB, batch by batch, and appends them to a results vector.
     *
     * Expected results - This test manually fetches sequences from a
     * SeqDB and computes the seeds for those sequences. The computed seeds
     * should be the same as the loaded ones.
     *
     * The GenerateMinimizers function is tested elsewhere.
     * The SeqDBReader::GetNext is also tested elsewhere.
    */

    // Input values.
    const std::string inSeqDB = PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/test-1.seqdb";
    const std::string inSeedDB =
        PacBio::PancakeTestsConfig::Data_Dir + "/seeddb-writer/test-1a.seeddb";
    int64_t batchSize = 10000;

    // Load the SeedDB.
    std::shared_ptr<PacBio::Pancake::SeedDBIndexCache> seedDBCache =
        PacBio::Pancake::LoadSeedDBIndexCache(inSeedDB);

    // Expected results - compute the minimizers from the sequences here.
    std::vector<PacBio::Pancake::FastaSequenceId> seqDBRecords = HelperLoadAllFromSeqDB(inSeqDB);
    std::vector<PacBio::Pancake::SequenceSeeds> computedSeeds = HelperComputeSequenceSeeds(
        seqDBRecords, seedDBCache->seedParams.KmerSize, seedDBCache->seedParams.MinimizerWindow,
        seedDBCache->seedParams.UseHPC, seedDBCache->seedParams.MaxHPCLen,
        seedDBCache->seedParams.UseRC);
    std::vector<std::vector<PacBio::Pancake::SequenceSeeds>> expected = {
        {computedSeeds[0], computedSeeds[1]},
        {computedSeeds[2]},
        {computedSeeds[3], computedSeeds[4]}};

    // Create a SeedDB reader.
    PacBio::Pancake::SeedDBReader reader(seedDBCache);

    // Collect all seeds from the SeedDB.
    std::vector<std::vector<PacBio::Pancake::SequenceSeeds>> results;
    std::vector<PacBio::Pancake::SequenceSeeds> batchRecords;
    while (reader.GetNextBatch(batchRecords, batchSize)) {
        results.emplace_back(batchRecords);
    }

    // Evaluate.
    EXPECT_EQ(expected, results);
}

TEST(SeedDBReader, GetNextBatch3)
{
    /*
     * This tests loading batches spread out through multiple files.
     * It loads all batches of 10kB, batch by batch, and appends them to a results vector.
     * Unlike the previous test, seeds for each sequence are stored in a separate
     * seeds file. Everything else is the same.
     *
     * Expected results - This test manually fetches sequences from a
     * SeqDB and computes the seeds for those sequences. The computed seeds
     * should be the same as the loaded ones.
     *
     * The GenerateMinimizers function is tested elsewhere.
     * The SeqDBReader::GetNext is also tested elsewhere.
    */

    // Input values.
    const std::string inSeqDB = PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/test-1.seqdb";
    const std::string inSeedDB =
        PacBio::PancakeTestsConfig::Data_Dir + "/seeddb-writer/test-1b.seeddb";
    int64_t batchSize = 10000;

    // Load the SeedDB.
    std::shared_ptr<PacBio::Pancake::SeedDBIndexCache> seedDBCache =
        PacBio::Pancake::LoadSeedDBIndexCache(inSeedDB);

    // Expected results - compute the minimizers from the sequences here.
    std::vector<PacBio::Pancake::FastaSequenceId> seqDBRecords = HelperLoadAllFromSeqDB(inSeqDB);
    std::vector<PacBio::Pancake::SequenceSeeds> computedSeeds = HelperComputeSequenceSeeds(
        seqDBRecords, seedDBCache->seedParams.KmerSize, seedDBCache->seedParams.MinimizerWindow,
        seedDBCache->seedParams.UseHPC, seedDBCache->seedParams.MaxHPCLen,
        seedDBCache->seedParams.UseRC);
    std::vector<std::vector<PacBio::Pancake::SequenceSeeds>> expected = {
        {computedSeeds[0], computedSeeds[1]},
        {computedSeeds[2]},
        {computedSeeds[3], computedSeeds[4]}};

    // Create a SeedDB reader.
    PacBio::Pancake::SeedDBReader reader(seedDBCache);

    // Collect all seeds from the SeedDB.
    std::vector<std::vector<PacBio::Pancake::SequenceSeeds>> results;
    std::vector<PacBio::Pancake::SequenceSeeds> batchRecords;
    while (reader.GetNextBatch(batchRecords, batchSize)) {
        results.emplace_back(batchRecords);
    }

    // Evaluate.
    EXPECT_EQ(expected, results);
}

TEST(SeedDBReader, GetBlock1)
{
    /*
     * Loads blocks in reverse, and accumulates the data into a vector.
     *
     * Expected results - This test manually fetches sequences from a
     * SeqDB and computes the seeds for those sequences. The computed seeds
     * should be the same as the loaded ones.
     *
     * The GenerateMinimizers function is tested elsewhere.
     * The SeqDBReader::GetNext is also tested elsewhere.
    */

    // Input values.
    const std::string inSeqDB = PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/test-1.seqdb";
    const std::string inSeedDB =
        PacBio::PancakeTestsConfig::Data_Dir + "/seeddb-writer/test-1a.seeddb";

    // Load the SeedDB.
    std::shared_ptr<PacBio::Pancake::SeedDBIndexCache> seedDBCache =
        PacBio::Pancake::LoadSeedDBIndexCache(inSeedDB);

    // Expected results - compute the minimizers from the sequences here.
    std::vector<PacBio::Pancake::FastaSequenceId> seqDBRecords = HelperLoadAllFromSeqDB(inSeqDB);
    std::vector<PacBio::Pancake::SequenceSeeds> computedSeeds = HelperComputeSequenceSeeds(
        seqDBRecords, seedDBCache->seedParams.KmerSize, seedDBCache->seedParams.MinimizerWindow,
        seedDBCache->seedParams.UseHPC, seedDBCache->seedParams.MaxHPCLen,
        seedDBCache->seedParams.UseRC);
    std::vector<std::vector<PacBio::Pancake::SequenceSeeds>> expected = {{computedSeeds[4]},
                                                                         {computedSeeds[3]},
                                                                         {computedSeeds[2]},
                                                                         {computedSeeds[1]},
                                                                         {computedSeeds[0]}};

    // Create a SeedDB reader.
    PacBio::Pancake::SeedDBReader reader(seedDBCache);

    // Collect all seeds from the SeedDB.
    std::vector<std::vector<PacBio::Pancake::SequenceSeeds>> results;
    std::vector<PacBio::Pancake::SequenceSeeds> blockRecords;
    for (int32_t i = static_cast<int32_t>(computedSeeds.size()) - 1; i >= 0; --i) {
        reader.GetBlock(blockRecords, i);
        results.emplace_back(blockRecords);
    }

    // Evaluate.
    EXPECT_EQ(expected, results);
}

TEST(SeedDBReader, GetBlock2)
{
    /*
     * Same as before, but tests loading seeds from multiple .seeds files.
     *
     * Loads blocks in reverse, and accumulates the data into a vector.
     *
     * Expected results - This test manually fetches sequences from a
     * SeqDB and computes the seeds for those sequences. The computed seeds
     * should be the same as the loaded ones.
     *
     * The GenerateMinimizers function is tested elsewhere.
     * The SeqDBReader::GetNext is also tested elsewhere.
    */

    // Input values.
    const std::string inSeqDB = PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/test-1.seqdb";
    const std::string inSeedDB =
        PacBio::PancakeTestsConfig::Data_Dir + "/seeddb-writer/test-1b.seeddb";

    // Load the SeedDB.
    std::shared_ptr<PacBio::Pancake::SeedDBIndexCache> seedDBCache =
        PacBio::Pancake::LoadSeedDBIndexCache(inSeedDB);

    // Expected results - compute the minimizers from the sequences here.
    std::vector<PacBio::Pancake::FastaSequenceId> seqDBRecords = HelperLoadAllFromSeqDB(inSeqDB);
    std::vector<PacBio::Pancake::SequenceSeeds> computedSeeds = HelperComputeSequenceSeeds(
        seqDBRecords, seedDBCache->seedParams.KmerSize, seedDBCache->seedParams.MinimizerWindow,
        seedDBCache->seedParams.UseHPC, seedDBCache->seedParams.MaxHPCLen,
        seedDBCache->seedParams.UseRC);
    std::vector<std::vector<PacBio::Pancake::SequenceSeeds>> expected = {{computedSeeds[4]},
                                                                         {computedSeeds[3]},
                                                                         {computedSeeds[2]},
                                                                         {computedSeeds[1]},
                                                                         {computedSeeds[0]}};

    // Create a SeedDB reader.
    PacBio::Pancake::SeedDBReader reader(seedDBCache);

    // Collect all seeds from the SeedDB.
    std::vector<std::vector<PacBio::Pancake::SequenceSeeds>> results;
    std::vector<PacBio::Pancake::SequenceSeeds> blockRecords;
    for (int32_t i = static_cast<int32_t>(computedSeeds.size()) - 1; i >= 0; --i) {
        reader.GetBlock(blockRecords, i);
        results.emplace_back(blockRecords);
    }

    // Evaluate.
    EXPECT_EQ(expected, results);
}

TEST(SeedDBReader, GetBlock3)
{
    /*
     * Tries to load a block out of bounds.
     * This should throw.
     *
     * The GenerateMinimizers function is tested elsewhere.
     * The SeqDBReader::GetNext is also tested elsewhere.
    */

    // Input values.
    const std::string inSeedDB =
        PacBio::PancakeTestsConfig::Data_Dir + "/seeddb-writer/test-1b.seeddb";

    // Load the SeedDB.
    std::shared_ptr<PacBio::Pancake::SeedDBIndexCache> seedDBCache =
        PacBio::Pancake::LoadSeedDBIndexCache(inSeedDB);

    // Create a SeedDB reader.
    PacBio::Pancake::SeedDBReader reader(seedDBCache);

    // Evaluate.
    EXPECT_THROW(
        {
            std::vector<PacBio::Pancake::SequenceSeeds> records;
            reader.GetBlock(records, -1);
        },
        std::runtime_error);

    EXPECT_THROW(
        {
            std::vector<PacBio::Pancake::SequenceSeeds> records;
            reader.GetBlock(records, 123);
        },
        std::runtime_error);
}

TEST(SeedDBReader, JumpTo1)
{
    /*
     * Uses the JumpTo function to position within the SeedDB, and then
     * tests that the location is correct by loading the next sequence
     * with GetNext.
     * This tests setting the JumpTo by sequence name.
     *
     * Expected results - This test manually fetches sequences from a
     * SeqDB and computes the seeds for those sequences. The computed seeds
     * should be the same as the loaded ones.
     *
     * The GenerateMinimizers function is tested elsewhere.
     * The SeqDBReader::GetNext is also tested elsewhere.
    */

    // Input values.
    const std::string inSeqDB = PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/test-1.seqdb";
    const std::string inSeedDB =
        PacBio::PancakeTestsConfig::Data_Dir + "/seeddb-writer/test-1a.seeddb";

    // Load the SeedDB.
    std::shared_ptr<PacBio::Pancake::SeedDBIndexCache> seedDBCache =
        PacBio::Pancake::LoadSeedDBIndexCache(inSeedDB);

    // Expected results - compute the minimizers from the sequences here.
    std::vector<PacBio::Pancake::FastaSequenceId> seqDBRecords = HelperLoadAllFromSeqDB(inSeqDB);
    std::vector<PacBio::Pancake::SequenceSeeds> expected = HelperComputeSequenceSeeds(
        seqDBRecords, seedDBCache->seedParams.KmerSize, seedDBCache->seedParams.MinimizerWindow,
        seedDBCache->seedParams.UseHPC, seedDBCache->seedParams.MaxHPCLen,
        seedDBCache->seedParams.UseRC);

    // Create a SeedDB reader.
    PacBio::Pancake::SeedDBReader reader(seedDBCache);

    // Reverse access all SequenceSeeds and fetch the ID from the SeedDB.
    // Check if they're the same.
    for (int32_t i = static_cast<int32_t>(expected.size()) - 1; i >= 0; --i) {
        const auto& expRecord = expected[i];

        reader.JumpTo(expRecord.Name());

        PacBio::Pancake::SequenceSeeds record;
        reader.GetNext(record);

        EXPECT_EQ(expRecord, record);
    }
}

TEST(SeedDBReader, JumpTo2)
{
    /*
     * Uses the JumpTo function to position within the SeedDB, and then
     * tests that the location is correct by loading the next sequence
     * with GetNext.
     * This tests setting the JumpTo by sequence ID.
     *
     * Expected results - This test manually fetches sequences from a
     * SeqDB and computes the seeds for those sequences. The computed seeds
     * should be the same as the loaded ones.
     *
     * The GenerateMinimizers function is tested elsewhere.
     * The SeqDBReader::GetNext is also tested elsewhere.
    */

    // Input values.
    const std::string inSeqDB = PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/test-1.seqdb";
    const std::string inSeedDB =
        PacBio::PancakeTestsConfig::Data_Dir + "/seeddb-writer/test-1a.seeddb";

    // Load the SeedDB.
    std::shared_ptr<PacBio::Pancake::SeedDBIndexCache> seedDBCache =
        PacBio::Pancake::LoadSeedDBIndexCache(inSeedDB);

    // Expected results - compute the minimizers from the sequences here.
    std::vector<PacBio::Pancake::FastaSequenceId> seqDBRecords = HelperLoadAllFromSeqDB(inSeqDB);
    std::vector<PacBio::Pancake::SequenceSeeds> expected = HelperComputeSequenceSeeds(
        seqDBRecords, seedDBCache->seedParams.KmerSize, seedDBCache->seedParams.MinimizerWindow,
        seedDBCache->seedParams.UseHPC, seedDBCache->seedParams.MaxHPCLen,
        seedDBCache->seedParams.UseRC);

    // Create a SeedDB reader.
    PacBio::Pancake::SeedDBReader reader(seedDBCache);

    // Reverse access all SequenceSeeds and fetch the ID from the SeedDB.
    // Check if they're the same.
    for (int32_t i = static_cast<int32_t>(expected.size()) - 1; i >= 0; --i) {
        const auto& expRecord = expected[i];

        reader.JumpTo(expRecord.Id());

        PacBio::Pancake::SequenceSeeds record;
        reader.GetNext(record);

        EXPECT_EQ(expRecord, record);
    }
}

TEST(SeedDBReader, JumpTo3)
{
    /*
     * Tries to jump to a sequence out of bounds.
     * This should throw.
     *
     * The GenerateMinimizers function is tested elsewhere.
     * The SeqDBReader::GetNext is also tested elsewhere.
    */

    // Input values.
    const std::string inSeedDB =
        PacBio::PancakeTestsConfig::Data_Dir + "/seeddb-writer/test-1b.seeddb";

    // Load the SeedDB.
    std::shared_ptr<PacBio::Pancake::SeedDBIndexCache> seedDBCache =
        PacBio::Pancake::LoadSeedDBIndexCache(inSeedDB);

    // Create a SeedDB reader.
    PacBio::Pancake::SeedDBReader reader(seedDBCache);

    // Evaluate.
    EXPECT_THROW(
        {
            std::vector<PacBio::Pancake::SequenceSeeds> records;
            reader.JumpTo(-1);
        },
        std::runtime_error);
    EXPECT_THROW(
        {
            std::vector<PacBio::Pancake::SequenceSeeds> records;
            reader.JumpTo(123);
        },
        std::runtime_error);
    EXPECT_THROW(
        {
            std::vector<PacBio::Pancake::SequenceSeeds> records;
            reader.JumpTo("nonexistent-name");
        },
        std::runtime_error);
}

TEST(SeedDBReader, MalformedSeedDB1)
{
    /*
     * The number of bytes for each sequence is wrong in the SeedDB, and does not match
     * the actual number of bytes.
     * The number of seeds for each sequence is correct though.
     * This should throw upon read.
     *
     * The GenerateMinimizers function is tested elsewhere.
     * The SeqDBReader::GetNext is also tested elsewhere.
    */

    // Input values.
    const std::string inSeedDB =
        PacBio::PancakeTestsConfig::Data_Dir + "/seeddb-writer/test-2a.seeddb";

    // Load the SeedDB.
    std::shared_ptr<PacBio::Pancake::SeedDBIndexCache> seedDBCache =
        PacBio::Pancake::LoadSeedDBIndexCache(inSeedDB);

    // Create a SeedDB reader.
    PacBio::Pancake::SeedDBReader reader(seedDBCache);

    // Test loading from the malformed SeedDB.
    EXPECT_THROW(
        {
            PacBio::Pancake::SequenceSeeds record;
            while (reader.GetNext(record))
                ;
        },
        std::runtime_error);
}

TEST(SeedDBReader, MalformedSeedDB2)
{
    /*
     * The number of seeds for each sequence is wrong in the SeedDB, and does not match
     * the actual number of seeds.
     * The number of bytes for each sequence is correct though.
     * This should throw upon read.
     *
     * The GenerateMinimizers function is tested elsewhere.
     * The SeqDBReader::GetNext is also tested elsewhere.
    */

    // Input values.
    const std::string inSeedDB =
        PacBio::PancakeTestsConfig::Data_Dir + "/seeddb-writer/test-2b.seeddb";

    // Load the SeedDB.
    std::shared_ptr<PacBio::Pancake::SeedDBIndexCache> seedDBCache =
        PacBio::Pancake::LoadSeedDBIndexCache(inSeedDB);

    // Create a SeedDB reader.
    PacBio::Pancake::SeedDBReader reader(seedDBCache);

    // Test loading from the malformed SeedDB.
    EXPECT_THROW(
        {
            PacBio::Pancake::SequenceSeeds record;
            while (reader.GetNext(record))
                ;
        },
        std::runtime_error);
}
