// Authors: Ivan Sovic

#include <PancakeTestData.h>
#include <gtest/gtest.h>
#include <pacbio/seqdb/SeqDBReader.h>
#include <pbbam/FastaReader.h>
#include <pbbam/FastaSequence.h>
#include <sstream>

std::vector<PacBio::BAM::FastaSequence> HelperLoadFasta(const std::string& inFasta)
{
    std::vector<PacBio::BAM::FastaSequence> ret;
    PacBio::BAM::FastaReader inReader{inFasta};
    PacBio::BAM::FastaSequence record;
    while (inReader.GetNext(record))
        ret.emplace_back(record);
    return ret;
}

std::string HelperLoadFastaAsString(const std::string& inFasta)
{
    std::ostringstream oss;
    auto records = HelperLoadFasta(inFasta);
    for (const auto& record : records)
        oss << ">" << record.Name() << "\n" << record.Bases() << "\n";
    return oss.str();
}

TEST(SeqDBReaderCompressed, GetNext)
{
    // Input values;
    const std::string inSeqFasta = PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/in.fasta";
    const std::string inSeqDB = PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/test-1.seqdb";

    // Load the expected sequences.
    // const std::string expected = HelperLoadFastaAsString(inSeqFasta);
    const auto expected = HelperLoadFasta(inSeqFasta);

    // Load the DB.
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBCache =
        PacBio::Pancake::LoadSeqDBIndexCache(inSeqDB);

    // Create a reader.
    PacBio::Pancake::SeqDBReader reader(seqDBCache);

    // Collect all sequences as a single string for comparison.
    std::vector<PacBio::BAM::FastaSequence> result;
    PacBio::BAM::FastaSequence record;
    while (reader.GetNext(record)) {
        result.emplace_back(record);
    }

    EXPECT_EQ(expected, result);
}

TEST(SeqDBReaderCompressed, GetSequenceByID1)
{
    /*
     * This tests accessing the sequences by their sequence ID.
    */

    // Input values;
    const std::string inSeqFasta = PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/in.fasta";
    const std::string inSeqDB = PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/test-1.seqdb";

    // Load the expected sequences.
    const auto expected = HelperLoadFasta(inSeqFasta);

    // Load the DB.
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBCache =
        PacBio::Pancake::LoadSeqDBIndexCache(inSeqDB);

    // Create a reader.
    PacBio::Pancake::SeqDBReader reader(seqDBCache);

    // Reverse access all sequences by ID and check if they're the same.
    for (int32_t i = static_cast<int32_t>(expected.size()) - 1; i >= 0; --i) {
        const auto& expRecord = expected[i];
        PacBio::BAM::FastaSequence record;
        reader.GetSequence(record, i);
        EXPECT_EQ(expRecord, record);
    }
}

TEST(SeqDBReaderCompressed, GetSequenceByID2)
{
    /*
     * Unlike the previous test (GetSequenceByID1) this loads a DB which as
     * permuted seqId values compared to their ordinal IDs. That means,
     * sequence with seqId = 3 appears as the first sequence in the file,
     * record with seqId = 0 is the second one, record with seqId = 4 is the third
     * one, and so on.
     * This tests random access via the sequence ID, which is not the same as the
     * ordinal ID.
    */

    // Input values;
    const std::string inSeqFasta = PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/in.fasta";
    const std::string inSeqDB = PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/test-2.seqdb";
    std::vector<int32_t> testSeqIds = {3, 0, 4, 2, 1};

    // Load the expected sequences.
    const auto fastaSeqs = HelperLoadFasta(inSeqFasta);
    std::vector<PacBio::BAM::FastaSequence> expected = {fastaSeqs[3], fastaSeqs[0], fastaSeqs[4],
                                                        fastaSeqs[2], fastaSeqs[1]};

    // Load the DB.
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBCache =
        PacBio::Pancake::LoadSeqDBIndexCache(inSeqDB);

    // Create a reader.
    PacBio::Pancake::SeqDBReader reader(seqDBCache);

    // Run the test.
    std::vector<PacBio::BAM::FastaSequence> results;
    for (const auto& seqId : testSeqIds) {
        PacBio::BAM::FastaSequence record;
        reader.GetSequence(record, seqId);
        results.emplace_back(record);
    }
    EXPECT_EQ(expected, results);
}

TEST(SeqDBReaderCompressed, GetSequenceByName)
{
    // Input values;
    const std::string inSeqFasta = PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/in.fasta";
    const std::string inSeqDB = PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/test-1.seqdb";

    // Load the expected sequences.
    const auto expected = HelperLoadFasta(inSeqFasta);

    // Load the DB.
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBCache =
        PacBio::Pancake::LoadSeqDBIndexCache(inSeqDB);

    // Create a reader.
    PacBio::Pancake::SeqDBReader reader(seqDBCache);

    // Reverse access all sequences by name and check if they're the same.
    for (int32_t i = static_cast<int32_t>(expected.size()) - 1; i >= 0; --i) {
        const auto& expRecord = expected[i];
        PacBio::BAM::FastaSequence record;
        reader.GetSequence(record, expRecord.Name());
        EXPECT_EQ(expRecord, record);
    }
}

TEST(SeqDBReaderCompressed, GetSequenceByIdNotExisting)
{
    // Input values;
    const std::string inSeqFasta = PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/in.fasta";
    const std::string inSeqDB = PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/test-1.seqdb";

    // Load the DB.
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBCache =
        PacBio::Pancake::LoadSeqDBIndexCache(inSeqDB);

    // Create a reader.
    PacBio::Pancake::SeqDBReader reader(seqDBCache);

    // Fetchin non-existent sequence names throws.
    {
        PacBio::BAM::FastaSequence record;
        try {
            reader.GetSequence(record, 10);
            ASSERT_TRUE(false);
        } catch (const std::runtime_error& e) {
        }
    }
    {
        PacBio::BAM::FastaSequence record;
        try {
            reader.GetSequence(record, -1);
            ASSERT_TRUE(false);
        } catch (const std::runtime_error& e) {
        }
    }
}

TEST(SeqDBReaderCompressed, GetSequenceByNameNotExisting)
{
    // Input values;
    const std::string inSeqFasta = PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/in.fasta";
    const std::string inSeqDB = PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/test-1.seqdb";

    // Load the DB.
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBCache =
        PacBio::Pancake::LoadSeqDBIndexCache(inSeqDB);

    // Create a reader.
    PacBio::Pancake::SeqDBReader reader(seqDBCache);

    // Fetching non-existent sequence names does no throw but returns false.
    // This is because it's hard for a user to know if the sequence is existing
    // or not, so hard to prevent throws.
    {
        PacBio::BAM::FastaSequence record;
        bool rv = reader.GetSequence(record, "some_nonexistent_name");
        ASSERT_FALSE(rv);
    }
}

TEST(SeqDBReaderCompressed, GetNextBatchAllAtOnce)
{
    // Input values;
    const int32_t batchSize = -1;
    const std::string inSeqFasta = PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/in.fasta";
    const std::string inSeqDB = PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/test-1.seqdb";

    // Load the expected sequences.
    const auto fastaSeqs = HelperLoadFasta(inSeqFasta);
    // Expect all sequences in the same block.
    const std::vector<std::vector<PacBio::BAM::FastaSequence>> expected = {
        {fastaSeqs[0], fastaSeqs[1], fastaSeqs[2], fastaSeqs[3], fastaSeqs[4]}};

    // Load the DB.
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBCache =
        PacBio::Pancake::LoadSeqDBIndexCache(inSeqDB);

    // Create a reader.
    PacBio::Pancake::SeqDBReader reader(seqDBCache);

    // Collect all sequences as a single string for comparison.
    std::vector<PacBio::BAM::FastaSequence> records;
    std::vector<std::vector<PacBio::BAM::FastaSequence>> results;
    while (reader.GetNextBatch(records, batchSize)) {
        results.emplace_back(records);
    }
    EXPECT_EQ(expected, results);
}

TEST(SeqDBReaderCompressed, GetNextBatchEverySeqSeparately)
{
    // Input values;
    const int32_t batchSize = 0;
    const std::string inSeqFasta = PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/in.fasta";
    const std::string inSeqDB = PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/test-1.seqdb";

    // Load the expected sequences.
    const auto fastaSeqs = HelperLoadFasta(inSeqFasta);
    // Expect all sequences in the same block.
    const std::vector<std::vector<PacBio::BAM::FastaSequence>> expected = {
        {fastaSeqs[0]}, {fastaSeqs[1]}, {fastaSeqs[2]}, {fastaSeqs[3]}, {fastaSeqs[4]}};

    // Load the DB.
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBCache =
        PacBio::Pancake::LoadSeqDBIndexCache(inSeqDB);

    // Create a reader.
    PacBio::Pancake::SeqDBReader reader(seqDBCache);

    // Collect all sequences as a single string for comparison.
    std::vector<PacBio::BAM::FastaSequence> records;
    std::vector<std::vector<PacBio::BAM::FastaSequence>> results;
    while (reader.GetNextBatch(records, batchSize)) {
        results.emplace_back(records);
    }
    EXPECT_EQ(expected, results);
}

TEST(SeqDBReaderCompressed, GetNextBatch10kbp)
{
    // Input values;
    const int32_t batchSize = 10000;
    const std::string inSeqFasta = PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/in.fasta";
    const std::string inSeqDB = PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/test-1.seqdb";

    // Load the expected sequences.
    const auto fastaSeqs = HelperLoadFasta(inSeqFasta);
    // Expect 3 blocks of sequences.
    const std::vector<std::vector<PacBio::BAM::FastaSequence>> expected = {
        {fastaSeqs[0], fastaSeqs[1]}, {fastaSeqs[2]}, {fastaSeqs[3], fastaSeqs[4]}};

    // Load the DB.
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBCache =
        PacBio::Pancake::LoadSeqDBIndexCache(inSeqDB);

    // Create a reader.
    PacBio::Pancake::SeqDBReader reader(seqDBCache);

    // Load sequences in batches of ~10kbp.
    std::vector<PacBio::BAM::FastaSequence> records;
    std::vector<std::vector<PacBio::BAM::FastaSequence>> results;
    while (reader.GetNextBatch(records, batchSize)) {
        results.emplace_back(records);
    }

    EXPECT_EQ(expected, results);
}

TEST(SeqDBReaderCompressed, GetBlockAllSeqsOneBlock)
{
    // Input values;
    const int32_t blockSize = -1;  // All sequences are in the same block.
    const std::string inSeqFasta = PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/in.fasta";
    const std::string inSeqDB = PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/test-1.seqdb";
    const std::vector<int32_t> testBlocks = {0};

    // Load the expected sequences.
    const auto fastaSeqs = HelperLoadFasta(inSeqFasta);
    // Expect 3 blocks of sequences.
    const std::vector<std::vector<PacBio::BAM::FastaSequence>> expected = {
        {fastaSeqs[0], fastaSeqs[1], fastaSeqs[2], fastaSeqs[3], fastaSeqs[4]}};

    // Load the DB.
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBCache =
        PacBio::Pancake::LoadSeqDBIndexCache(inSeqDB, blockSize);

    // Create a reader.
    PacBio::Pancake::SeqDBReader reader(seqDBCache);

    // Load sequences in batches of ~10kbp.
    std::vector<PacBio::BAM::FastaSequence> records;
    std::vector<std::vector<PacBio::BAM::FastaSequence>> results;
    for (const auto& blockId : testBlocks) {
        reader.GetBlock(records, blockId);
        results.emplace_back(records);
    }

    EXPECT_EQ(expected, results);
}

TEST(SeqDBReaderCompressed, GetBlockEachSeqSeparateBlock)
{
    // Input values;
    const int32_t blockSize = 0;  // All sequences are in the same block.
    const std::string inSeqFasta = PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/in.fasta";
    const std::string inSeqDB = PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/test-1.seqdb";
    const std::vector<int32_t> testBlocks = {0, 1, 2, 3, 4};

    // Load the expected sequences.
    const auto fastaSeqs = HelperLoadFasta(inSeqFasta);
    // Expect 3 blocks of sequences.
    const std::vector<std::vector<PacBio::BAM::FastaSequence>> expected = {
        {fastaSeqs[0]}, {fastaSeqs[1]}, {fastaSeqs[2]}, {fastaSeqs[3]}, {fastaSeqs[4]}};

    // Load the DB.
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBCache =
        PacBio::Pancake::LoadSeqDBIndexCache(inSeqDB, blockSize);

    // Create a reader.
    PacBio::Pancake::SeqDBReader reader(seqDBCache);

    // Load sequences in batches of ~10kbp.
    std::vector<PacBio::BAM::FastaSequence> records;
    std::vector<std::vector<PacBio::BAM::FastaSequence>> results;
    for (const auto& blockId : testBlocks) {
        reader.GetBlock(records, blockId);
        results.emplace_back(records);
    }

    EXPECT_EQ(expected, results);
}

TEST(SeqDBReaderCompressed, GetBlock10kbp)
{
    // Input values;
    const int32_t blockSize = 10000;  // 10kbp of compressed sequences.
    const std::string inSeqFasta = PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/in.fasta";
    const std::string inSeqDB = PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/test-1.seqdb";
    const std::vector<int32_t> testBlocks = {0, 1};

    // Load the expected sequences.
    const auto fastaSeqs = HelperLoadFasta(inSeqFasta);
    // Expect 3 blocks of sequences.
    const std::vector<std::vector<PacBio::BAM::FastaSequence>> expected = {
        {fastaSeqs[0], fastaSeqs[1], fastaSeqs[2]}, {fastaSeqs[3], fastaSeqs[4]}};

    // Load the DB.
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBCache =
        PacBio::Pancake::LoadSeqDBIndexCache(inSeqDB, blockSize);

    // Create a reader.
    PacBio::Pancake::SeqDBReader reader(seqDBCache);

    // Load sequences in batches of ~10kbp.
    std::vector<PacBio::BAM::FastaSequence> records;
    std::vector<std::vector<PacBio::BAM::FastaSequence>> results;
    for (const auto& blockId : testBlocks) {
        reader.GetBlock(records, blockId);
        results.emplace_back(records);
    }

    EXPECT_EQ(expected, results);
}

TEST(SeqDBReaderCompressed, GetBlockOutOfBounds)
{
    // Input values;
    int32_t blockSize = 10000;  // 10kbp of compressed sequences.
    const std::string inSeqFasta = PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/in.fasta";
    const std::string inSeqDB = PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/test-1.seqdb";
    const std::vector<int32_t> testBlocks = {-1, 10};

    // Load the DB.
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBCache =
        PacBio::Pancake::LoadSeqDBIndexCache(inSeqDB, blockSize);

    // Create a reader.
    PacBio::Pancake::SeqDBReader reader(seqDBCache);

    // Load sequences in batches of ~10kbp.
    std::vector<PacBio::BAM::FastaSequence> records;
    for (const auto& blockId : testBlocks) {
        try {
            reader.GetBlock(records, blockId);
            ASSERT_TRUE(false);
        } catch (const std::runtime_error& e) {
        }
    }
}

TEST(SeqDBReaderCompressed, JumpToByName)
{
    // Input values;
    const std::string inSeqFasta = PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/in.fasta";
    const std::string inSeqDB = PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/test-1.seqdb";
    const std::vector<std::string> testNames = {
        "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3981/0_5105",   // 3
        "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3820/0_24292",  // 2
        "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/4028/0_19001",  // 4
        "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3005/0_5852",   // 0
        "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3414/0_11983",  // 1
    };

    // Load the expected sequences.
    const auto fastaSeqs = HelperLoadFasta(inSeqFasta);
    const std::vector<PacBio::BAM::FastaSequence> expected = {
        fastaSeqs[3], fastaSeqs[2], fastaSeqs[4], fastaSeqs[0], fastaSeqs[1]};

    // Load the DB.
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBCache =
        PacBio::Pancake::LoadSeqDBIndexCache(inSeqDB);

    // Create a reader.
    PacBio::Pancake::SeqDBReader reader(seqDBCache);

    // Set the positon to the next sequence in the order specified by testNames.
    // The GetNext should then load that particular sequence.
    std::vector<PacBio::BAM::FastaSequence> results;
    PacBio::BAM::FastaSequence record;
    for (const auto& name : testNames) {
        reader.JumpTo(name);
        reader.GetNext(record);
        results.emplace_back(record);
    }

    EXPECT_EQ(expected, results);
}

TEST(SeqDBReaderCompressed, JumpToById)
{
    // Input values;
    const std::string inSeqFasta = PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/in.fasta";
    const std::string inSeqDB = PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/test-1.seqdb";
    const std::vector<int32_t> testIds = {3, 2, 4, 0, 1};

    // Load the expected sequences.
    const auto fastaSeqs = HelperLoadFasta(inSeqFasta);
    const std::vector<PacBio::BAM::FastaSequence> expected = {
        fastaSeqs[3], fastaSeqs[2], fastaSeqs[4], fastaSeqs[0], fastaSeqs[1]};

    // Load the DB.
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBCache =
        PacBio::Pancake::LoadSeqDBIndexCache(inSeqDB);

    // Create a reader.
    PacBio::Pancake::SeqDBReader reader(seqDBCache);

    // Set the positon to the next sequence in the order specified by testNames.
    // The GetNext should then load that particular sequence.
    std::vector<PacBio::BAM::FastaSequence> results;
    PacBio::BAM::FastaSequence record;
    for (const auto& id : testIds) {
        reader.JumpTo(id);
        reader.GetNext(record);
        results.emplace_back(record);
    }

    EXPECT_EQ(expected, results);
}

TEST(SeqDBReaderUncompressed, GetNext)
{
    /*
     * This is an equivalent test to TEST(SeqDBReaderCompressed, GetNext), but the
     * test-3.seqdb contains uncompressed sequences.
     * This is enough to prove that uncompressed loading works, because the rest
     * of the functionality is the same in both cases.
    */

    // Input values;
    const std::string inSeqFasta = PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/in.fasta";
    const std::string inSeqDB = PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/test-3.seqdb";

    // Load the expected sequences.
    // const std::string expected = HelperLoadFastaAsString(inSeqFasta);
    const auto expected = HelperLoadFasta(inSeqFasta);

    // Load the DB.
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBCache =
        PacBio::Pancake::LoadSeqDBIndexCache(inSeqDB);

    // Create a reader.
    PacBio::Pancake::SeqDBReader reader(seqDBCache);

    // Collect all sequences as a single string for comparison.
    std::vector<PacBio::BAM::FastaSequence> result;
    PacBio::BAM::FastaSequence record;
    while (reader.GetNext(record)) {
        result.emplace_back(record);
    }

    EXPECT_EQ(expected, result);
}

TEST(SeqDBReader, MalformedBytes)
{
    /*
     * The sequence length in these DBs is wrong - the DBs specify 10000 bp long
     * sequences, but in reality, the sequence is of length 5852 bases.
     * The test-4.seqdb runs the parser with the decompression step, and
     * test-5.seqdb without it.
     * Both DBs reuse existing data from previous tests (test-1 for compressed, and
     * test-5 for uncompressed).
    */

    // Input values;
    const std::string inSeqFasta = PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/in.fasta";
    const std::vector<std::string> inSeqDBList = {
        PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/test-4.seqdb",
        PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/test-5.seqdb"};

    // Run tests.
    for (const auto& inSeqDB : inSeqDBList) {
        // Load the DB.
        std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBCache =
            PacBio::Pancake::LoadSeqDBIndexCache(inSeqDB);

        // Create a reader.
        PacBio::Pancake::SeqDBReader reader(seqDBCache);

        // Collect all sequences as a single string for comparison.
        PacBio::BAM::FastaSequence record;
        while (true) {
            bool rv = false;
            try {
                rv = reader.GetNext(record);
                ASSERT_TRUE(false);
            } catch (const std::runtime_error& e) {
            }
            if (!rv) break;
        }
    }
}
