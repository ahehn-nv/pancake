#include <gtest/gtest.h>
#include <pacbio/seeddb/Minimizers.h>
#include <pacbio/seeddb/Seed.h>

void HelperTestGenerateMinimizers(const std::string& seq, int32_t seqId, int32_t k, int32_t w,
                                  int32_t space, bool useHPC, int32_t maxHPCLen, bool useRC,
                                  int32_t expectedRv, const std::vector<__int128>& expectedSeeds)
{
    // Run the unit under test.
    const uint8_t* seqData = reinterpret_cast<const uint8_t*>(seq.data());
    int32_t seqLen = seq.size();
    std::vector<__int128> seeds;
    int rv = PacBio::Pancake::SeedDB::GenerateMinimizers(seeds, seqData, seqLen, 0, seqId, k, w,
                                                         space, useRC, useHPC, maxHPCLen);

    // Compare the results.
    EXPECT_EQ(expectedRv, rv);
    EXPECT_EQ(expectedSeeds, seeds);
}

TEST(GenerateMinimizers, SmallTest1)
{
    // Inputs.
    const std::string seq = "AAAAAAAAAA";
    const int32_t seqId = 123;
    const int32_t k = 5;
    const int32_t w = 1;
    const int32_t space = 0;
    const bool useHPC = false;
    const int32_t maxHPCLen = 10;
    const bool useRC = true;

    // Expected results.
    const int32_t expectedRv = 0;
    const std::vector<__int128> expectedSeeds = {
        PacBio::Pancake::SeedDB::Seed::Encode(0, seqId, 0, false),
        PacBio::Pancake::SeedDB::Seed::Encode(0, seqId, 1, false),
        PacBio::Pancake::SeedDB::Seed::Encode(0, seqId, 2, false),
        PacBio::Pancake::SeedDB::Seed::Encode(0, seqId, 3, false),
        PacBio::Pancake::SeedDB::Seed::Encode(0, seqId, 4, false),
        PacBio::Pancake::SeedDB::Seed::Encode(0, seqId, 5, false),
    };

    HelperTestGenerateMinimizers(seq, seqId, k, w, space, useHPC, maxHPCLen, useRC, expectedRv,
                                 expectedSeeds);
}

TEST(GenerateMinimizers, SmallTest2)
{
    // Inputs.
    const std::string seq = "CAAAAAAAAA";
    const int32_t seqId = 123;
    const int32_t k = 5;
    const int32_t w = 1;
    const int32_t space = 0;
    const bool useHPC = false;
    const int32_t maxHPCLen = 10;
    const bool useRC = true;

    // Expected results.
    const int32_t expectedRv = 0;
    const std::vector<__int128> expectedSeeds = {
        PacBio::Pancake::SeedDB::Seed::Encode(256, seqId, 0, false),
        PacBio::Pancake::SeedDB::Seed::Encode(0, seqId, 1, false),
        PacBio::Pancake::SeedDB::Seed::Encode(0, seqId, 2, false),
        PacBio::Pancake::SeedDB::Seed::Encode(0, seqId, 3, false),
        PacBio::Pancake::SeedDB::Seed::Encode(0, seqId, 4, false),
        PacBio::Pancake::SeedDB::Seed::Encode(0, seqId, 5, false),
    };

    HelperTestGenerateMinimizers(seq, seqId, k, w, space, useHPC, maxHPCLen, useRC, expectedRv,
                                 expectedSeeds);
}

TEST(GenerateMinimizers, SmallTest3)
{
    // Inputs.
    const std::string seq = "AGCTTTTCATTCTGACTGCA";
    const int32_t seqId = 123;
    const int32_t k = 5;
    const int32_t w = 1;
    const int32_t space = 0;
    const bool useHPC = false;
    const int32_t maxHPCLen = 10;
    const bool useRC = true;

    // Expected results.
    const int32_t expectedRv = 0;
    const std::vector<__int128> expectedSeeds = {
        PacBio::Pancake::SeedDB::Seed::Encode(39, seqId, 0, true),
        PacBio::Pancake::SeedDB::Seed::Encode(9, seqId, 1, true),
        PacBio::Pancake::SeedDB::Seed::Encode(2, seqId, 2, true),
        PacBio::Pancake::SeedDB::Seed::Encode(512, seqId, 3, true),
        PacBio::Pancake::SeedDB::Seed::Encode(896, seqId, 4, true),
        PacBio::Pancake::SeedDB::Seed::Encode(224, seqId, 5, true),
        PacBio::Pancake::SeedDB::Seed::Encode(56, seqId, 6, true),
        PacBio::Pancake::SeedDB::Seed::Encode(317, seqId, 7, false),
        PacBio::Pancake::SeedDB::Seed::Encode(131, seqId, 8, true),
        PacBio::Pancake::SeedDB::Seed::Encode(288, seqId, 9, true),
        PacBio::Pancake::SeedDB::Seed::Encode(840, seqId, 10, true),
        PacBio::Pancake::SeedDB::Seed::Encode(481, seqId, 11, false),
        PacBio::Pancake::SeedDB::Seed::Encode(180, seqId, 12, true),
        PacBio::Pancake::SeedDB::Seed::Encode(301, seqId, 13, true),
        PacBio::Pancake::SeedDB::Seed::Encode(121, seqId, 14, false),
        PacBio::Pancake::SeedDB::Seed::Encode(484, seqId, 15, false),
    };

    HelperTestGenerateMinimizers(seq, seqId, k, w, space, useHPC, maxHPCLen, useRC, expectedRv,
                                 expectedSeeds);
}

TEST(GenerateMinimizers, SmallTest4)
{
    // Inputs.
    const std::string seq = "TCTGACTGCA";
    const int32_t seqId = 123;
    const int32_t k = 5;
    const int32_t w = 1;
    const int32_t space = 0;
    const bool useHPC = false;
    const int32_t maxHPCLen = 10;
    const bool useRC = true;

    // Expected results.
    const int32_t expectedRv = 0;
    const std::vector<__int128> expectedSeeds = {
        PacBio::Pancake::SeedDB::Seed::Encode(840, seqId, 0, true),
        PacBio::Pancake::SeedDB::Seed::Encode(481, seqId, 1, false),
        PacBio::Pancake::SeedDB::Seed::Encode(180, seqId, 2, true),
        PacBio::Pancake::SeedDB::Seed::Encode(301, seqId, 3, true),
        PacBio::Pancake::SeedDB::Seed::Encode(121, seqId, 4, false),
        PacBio::Pancake::SeedDB::Seed::Encode(484, seqId, 5, false),
    };

    HelperTestGenerateMinimizers(seq, seqId, k, w, space, useHPC, maxHPCLen, useRC, expectedRv,
                                 expectedSeeds);
}

TEST(GenerateMinimizers, SmallTest5)
{
    // Inputs.
    const std::string seq =
        "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC";
    const int32_t seqId = 123;
    const int32_t k = 15;
    const int32_t w = 5;
    const int32_t space = 0;
    const bool useHPC = false;
    const int32_t maxHPCLen = 10;
    const bool useRC = true;

    // Expected results.
    const int32_t expectedRv = 0;
    const std::vector<__int128> expectedSeeds = {
        PacBio::Pancake::SeedDB::Seed::Encode(167726968, seqId, 0, false),
        PacBio::Pancake::SeedDB::Seed::Encode(189282306, seqId, 2, true),
        PacBio::Pancake::SeedDB::Seed::Encode(315756032, seqId, 3, true),
        PacBio::Pancake::SeedDB::Seed::Encode(259553542, seqId, 8, false),
        PacBio::Pancake::SeedDB::Seed::Encode(364792648, seqId, 10, true),
        PacBio::Pancake::SeedDB::Seed::Encode(126904899, seqId, 14, false),
        PacBio::Pancake::SeedDB::Seed::Encode(80713150, seqId, 18, true),
        PacBio::Pancake::SeedDB::Seed::Encode(27856109, seqId, 19, false),
        PacBio::Pancake::SeedDB::Seed::Encode(76041465, seqId, 24, true),
        PacBio::Pancake::SeedDB::Seed::Encode(54228923, seqId, 26, false),
        PacBio::Pancake::SeedDB::Seed::Encode(55644705, seqId, 31, true),
        PacBio::Pancake::SeedDB::Seed::Encode(502196464, seqId, 33, false),
        PacBio::Pancake::SeedDB::Seed::Encode(518950656, seqId, 35, false),
        PacBio::Pancake::SeedDB::Seed::Encode(536855377, seqId, 39, true),
        PacBio::Pancake::SeedDB::Seed::Encode(503315509, seqId, 41, true),
        PacBio::Pancake::SeedDB::Seed::Encode(125828877, seqId, 42, true),
        PacBio::Pancake::SeedDB::Seed::Encode(74973168, seqId, 44, true),
        PacBio::Pancake::SeedDB::Seed::Encode(35767, seqId, 46, false),
        PacBio::Pancake::SeedDB::Seed::Encode(143070, seqId, 47, false),
        PacBio::Pancake::SeedDB::Seed::Encode(572280, seqId, 48, false),
        PacBio::Pancake::SeedDB::Seed::Encode(2289123, seqId, 49, false),
        PacBio::Pancake::SeedDB::Seed::Encode(9156492, seqId, 50, false),
        PacBio::Pancake::SeedDB::Seed::Encode(36625970, seqId, 51, false),
    };

    HelperTestGenerateMinimizers(seq, seqId, k, w, space, useHPC, maxHPCLen, useRC, expectedRv,
                                 expectedSeeds);
}

TEST(GenerateMinimizers, SmallTest6)
{
    // Inputs.
    const std::string seq = "AGCTTTTCATTCTGACTGCA";
    const int32_t seqId = 123;
    const int32_t k = 5;
    const int32_t w = 4;
    const int32_t space = 0;
    const bool useHPC = false;
    const int32_t maxHPCLen = 10;
    const bool useRC = true;

    // Expected results.
    const int32_t expectedRv = 0;
    const std::vector<__int128> expectedSeeds = {
        PacBio::Pancake::SeedDB::Seed::Encode(2, seqId, 2, true),
        PacBio::Pancake::SeedDB::Seed::Encode(56, seqId, 6, true),
        PacBio::Pancake::SeedDB::Seed::Encode(131, seqId, 8, true),
        PacBio::Pancake::SeedDB::Seed::Encode(180, seqId, 12, true),
        PacBio::Pancake::SeedDB::Seed::Encode(121, seqId, 14, false),
    };

    HelperTestGenerateMinimizers(seq, seqId, k, w, space, useHPC, maxHPCLen, useRC, expectedRv,
                                 expectedSeeds);
}

TEST(GenerateMinimizers, SmallTest7WithNBases)
{
    // Inputs.
    const std::string seq = "AGCTTTTCATTCTGACTGCANNNACTNNNNNAGCTTTTCATTCTGACTGCA";
    const int32_t seqId = 123;
    const int32_t k = 5;
    const int32_t w = 4;
    const int32_t space = 0;
    const bool useHPC = false;
    const int32_t maxHPCLen = 10;
    const bool useRC = true;

    // Expected results.
    const int32_t expectedRv = 0;
    const std::vector<__int128> expectedSeeds = {
        PacBio::Pancake::SeedDB::Seed::Encode(2, seqId, 2, true),
        PacBio::Pancake::SeedDB::Seed::Encode(56, seqId, 6, true),
        PacBio::Pancake::SeedDB::Seed::Encode(131, seqId, 8, true),
        PacBio::Pancake::SeedDB::Seed::Encode(180, seqId, 12, true),
        PacBio::Pancake::SeedDB::Seed::Encode(121, seqId, 14, false),
        PacBio::Pancake::SeedDB::Seed::Encode(2, seqId, 33, true),
        PacBio::Pancake::SeedDB::Seed::Encode(56, seqId, 37, true),
        PacBio::Pancake::SeedDB::Seed::Encode(131, seqId, 39, true),
        PacBio::Pancake::SeedDB::Seed::Encode(180, seqId, 43, true),
        PacBio::Pancake::SeedDB::Seed::Encode(121, seqId, 45, false),
    };

    HelperTestGenerateMinimizers(seq, seqId, k, w, space, useHPC, maxHPCLen, useRC, expectedRv,
                                 expectedSeeds);
}

TEST(GenerateMinimizers, SeedSize32BasePairs_PolyA)
{
    // Inputs.
    const std::string seq = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
    const int32_t seqId = 0;
    const int32_t k = 32;
    const int32_t w = 1;
    const bool useHPC = false;
    const int32_t maxHPCLen = 10;
    const bool useRC = true;

    // Expected results.
    const int32_t expectedRv = 0;
    const std::vector<__int128> expectedSeeds = {
        PacBio::Pancake::SeedDB::Seed::Encode(0, seqId, 0, false),
    };

    HelperTestGenerateMinimizers(seq, seqId, k, w, useHPC, maxHPCLen, useRC, expectedRv,
                                 expectedSeeds);
}

TEST(GenerateMinimizers, SeedSize32BasePairs_PolyT_WithRC)
{
    // Inputs.
    const std::string seq = "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";
    const int32_t seqId = 0;
    const int32_t k = 32;
    const int32_t w = 1;
    const bool useHPC = false;
    const int32_t maxHPCLen = 10;
    const bool useRC = true;

    // Expected results.
    const int32_t expectedRv = 0;
    const std::vector<__int128> expectedSeeds = {
        PacBio::Pancake::SeedDB::Seed::Encode(0, seqId, 0, true),
    };

    HelperTestGenerateMinimizers(seq, seqId, k, w, useHPC, maxHPCLen, useRC, expectedRv,
                                 expectedSeeds);
}

TEST(GenerateMinimizers, SeedSize32BasePairs_PolyT_OnlyFWD)
{
    // Inputs.
    const std::string seq = "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";
    const int32_t seqId = 0;
    const int32_t k = 32;
    const int32_t w = 1;
    const bool useHPC = false;
    const int32_t maxHPCLen = 10;
    const bool useRC = false;

    // Expected results.
    const int32_t expectedRv = 0;
    const std::vector<__int128> expectedSeeds = {
        PacBio::Pancake::SeedDB::Seed::Encode(0xFFFFFFFFFFFFFFFF, seqId, 0, false),
    };

    HelperTestGenerateMinimizers(seq, seqId, k, w, useHPC, maxHPCLen, useRC, expectedRv,
                                 expectedSeeds);
}

TEST(GenerateMinimizers, HPC1)
{
    // Inputs.
    const std::string seq = "ACGTTTTG";
    const int32_t seqId = 123;
    const int32_t k = 5;
    const int32_t w = 1;
    const int32_t space = 0;
    const bool useHPC = true;
    const int32_t maxHPCLen = 5;
    const bool useRC = true;

    // Expected results.
    const int32_t expectedRv = 0;
    const std::vector<__int128> expectedSeeds = {
        PacBio::Pancake::SeedDB::Seed::Encode(110, seqId, 0, false),
    };

    HelperTestGenerateMinimizers(seq, seqId, k, w, space, useHPC, maxHPCLen, useRC, expectedRv,
                                 expectedSeeds);
}

TEST(GenerateMinimizers, HPC2)
{
    // This test is a bit more comprehensive. It does not explicitly specify inputs and outputs, but performs
    // two runs with different settings and inputs, and expects the same results.
    // First run - use a sequence with homopolymer runs, and use the homopolymer compression options to generate the seeds.
    // Second run - use a sequence with manually compressed homopolymers, and compute minimizers without the compression option.
    // Both runs should result in identical sets of seeds, except for the positions (because one sequence is shorter than the other.

    // Inputs.
    const std::string seqWithHP =
        "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC";
    // Manually removed homopolymers from the sequence (manual HP compression).
    const std::string seqWithoutHP = "AGCTCATCTGACTGCACGCATATGTCTCTGTGTGATAGAGTGTCTGATAGCAGC";
    const int32_t seqId = 123;
    const int32_t k = 5;
    const int32_t w = 1;
    const int32_t space = 0;
    const int32_t maxHPCLen = 5;
    const bool useRC = true;

    // Run the unit under test WITH homopolymer compression for the sequence with homopolymers.
    std::vector<__int128> seedsWithHP;
    int32_t rvWithHP = 0;
    {
        const uint8_t* seqData = reinterpret_cast<const uint8_t*>(seqWithHP.data());
        int32_t seqLen = seqWithHP.size();
        rvWithHP = PacBio::Pancake::SeedDB::GenerateMinimizers(
            seedsWithHP, seqData, seqLen, 0, seqId, k, w, space, useRC, true, maxHPCLen);
    }

    // Run the unit under test WITHOUT homopolymer compression on the sequence with no homopolymers.
    std::vector<__int128> seedsNoHP;
    int32_t rvNoHP = 0;
    {
        const uint8_t* seqData = reinterpret_cast<const uint8_t*>(seqWithoutHP.data());
        int32_t seqLen = seqWithoutHP.size();
        rvNoHP = PacBio::Pancake::SeedDB::GenerateMinimizers(seedsNoHP, seqData, seqLen, 0, seqId,
                                                             k, w, space, useRC, false, maxHPCLen);
    }

    // Check the return values.
    EXPECT_EQ(rvWithHP, 0);
    EXPECT_EQ(rvNoHP, 0);

    // There should be the same amount of seeds in both vectors.
    EXPECT_EQ(seedsWithHP.size(), seedsNoHP.size());

    // Need to check each seed manually, because positions are different.
    // Both seed vectors should have the same seeds, just at different positions.
    for (size_t i = 0; i < seedsWithHP.size(); ++i) {
        auto seedWithHP = PacBio::Pancake::SeedDB::Seed(seedsWithHP[i]);
        auto seedNoHP = PacBio::Pancake::SeedDB::Seed(seedsNoHP[i]);
        EXPECT_EQ(seedWithHP.key, seedNoHP.key);
        EXPECT_EQ(seedWithHP.flag, seedNoHP.flag);
    }
}
