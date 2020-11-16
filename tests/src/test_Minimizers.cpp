#include <gtest/gtest.h>
#include <pacbio/pancake/Minimizers.h>
#include <pacbio/pancake/Seed.h>
#include <pacbio/util/CommonTypes.h>
// #include <iostream>

using namespace PacBio::Pancake;

void HelperTestGenerateMinimizers(const std::string& seq, int32_t seqId, int32_t k, int32_t w,
                                  int32_t space, bool useHPC, int32_t maxHPCLen, bool useRC,
                                  int32_t expectedRv,
                                  const std::vector<PacBio::Pancake::Int128t>& expectedSeeds,
                                  bool expectedThrow)
{
    // Run the unit under test.
    const uint8_t* seqData = reinterpret_cast<const uint8_t*>(seq.data());
    int32_t seqLen = seq.size();
    std::vector<PacBio::Pancake::Int128t> seeds;

    if (expectedThrow) {
        EXPECT_THROW(
            {
                PacBio::Pancake::SeedDB::GenerateMinimizers(seeds, seqData, seqLen, 0, seqId, k, w,
                                                            space, useRC, useHPC, maxHPCLen);
            },
            std::runtime_error);

    } else {
        const int rv = PacBio::Pancake::SeedDB::GenerateMinimizers(
            seeds, seqData, seqLen, 0, seqId, k, w, space, useRC, useHPC, maxHPCLen);

        // std::cerr << "Results:\n";
        // for (const auto& val : seeds) {
        //     auto s = PacBio::Pancake::SeedDB::Seed(val);
        //     // std::cerr << s.Verbose() << "\n";
        //     std::cerr << "PacBio::Pancake::SeedDB::Seed::Encode(" << s.key << ", " << s.span << ", "
        //               << s.seqID << ", " << s.pos << ", " << s.seqRev << "),\n";
        // }

        // std::cerr << "Expected:\n";
        // for (const auto& val : expectedSeeds) {
        //     auto s = PacBio::Pancake::SeedDB::Seed(val);
        //     // std::cerr << s.Verbose() << "\n";
        //     std::cerr << "PacBio::Pancake::SeedDB::Seed::Encode(" << s.key << ", " << s.span << ", "
        //               << s.seqID << ", " << s.pos << ", " << s.seqRev << "),\n";
        // }

        // Compare the results.
        EXPECT_EQ(expectedRv, rv);
        EXPECT_EQ(expectedSeeds, seeds);
    }
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

    // Helper function.
    const uint64_t mask = SeedDB::ComputeKmerMask(k);
    auto Hash = [&](uint64_t val) { return SeedDB::InvertibleHash(val, mask); };

    // Expected results.
    const int32_t expectedRv = 0;
    const bool expectedThrow = false;
    const std::vector<PacBio::Pancake::Int128t> expectedSeeds = {
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(0), k, seqId, 0, false),
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(0), k, seqId, 1, false),
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(0), k, seqId, 2, false),
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(0), k, seqId, 3, false),
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(0), k, seqId, 4, false),
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(0), k, seqId, 5, false),
    };

    HelperTestGenerateMinimizers(seq, seqId, k, w, space, useHPC, maxHPCLen, useRC, expectedRv,
                                 expectedSeeds, expectedThrow);
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

    // Helper function.
    const uint64_t mask = SeedDB::ComputeKmerMask(k);
    auto Hash = [&](uint64_t val) { return SeedDB::InvertibleHash(val, mask); };

    // Expected results.
    const int32_t expectedRv = 0;
    const bool expectedThrow = false;
    const std::vector<PacBio::Pancake::Int128t> expectedSeeds = {
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(256), k, seqId, 0, false),
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(0), k, seqId, 1, false),
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(0), k, seqId, 2, false),
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(0), k, seqId, 3, false),
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(0), k, seqId, 4, false),
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(0), k, seqId, 5, false),
    };

    HelperTestGenerateMinimizers(seq, seqId, k, w, space, useHPC, maxHPCLen, useRC, expectedRv,
                                 expectedSeeds, expectedThrow);
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

    // Helper function.
    const uint64_t mask = SeedDB::ComputeKmerMask(k);
    auto Hash = [&](uint64_t val) { return SeedDB::InvertibleHash(val, mask); };

    // Expected results.
    const int32_t expectedRv = 0;
    const bool expectedThrow = false;
    const std::vector<PacBio::Pancake::Int128t> expectedSeeds = {
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(39), k, seqId, 0, true),
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(9), k, seqId, 1, true),
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(2), k, seqId, 2, true),
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(512), k, seqId, 3, true),
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(896), k, seqId, 4, true),
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(224), k, seqId, 5, true),
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(56), k, seqId, 6, true),
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(317), k, seqId, 7, false),
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(131), k, seqId, 8, true),
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(288), k, seqId, 9, true),
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(840), k, seqId, 10, true),
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(481), k, seqId, 11, false),
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(180), k, seqId, 12, true),
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(301), k, seqId, 13, true),
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(121), k, seqId, 14, false),
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(484), k, seqId, 15, false),
    };

    HelperTestGenerateMinimizers(seq, seqId, k, w, space, useHPC, maxHPCLen, useRC, expectedRv,
                                 expectedSeeds, expectedThrow);
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

    // Helper function.
    const uint64_t mask = SeedDB::ComputeKmerMask(k);
    auto Hash = [&](uint64_t val) { return SeedDB::InvertibleHash(val, mask); };

    // Expected results.
    const int32_t expectedRv = 0;
    const bool expectedThrow = false;
    const std::vector<PacBio::Pancake::Int128t> expectedSeeds = {
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(840), k, seqId, 0, true),
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(481), k, seqId, 1, false),
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(180), k, seqId, 2, true),
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(301), k, seqId, 3, true),
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(121), k, seqId, 4, false),
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(484), k, seqId, 5, false),
    };

    HelperTestGenerateMinimizers(seq, seqId, k, w, space, useHPC, maxHPCLen, useRC, expectedRv,
                                 expectedSeeds, expectedThrow);
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
    const bool expectedThrow = false;
    /*
     * Without the invertible hash function, these would be the generated seeds.
    */
    // const std::vector<PacBio::Pancake::Int128t> expectedSeeds = {
    //     PacBio::Pancake::SeedDB::Seed::Encode(Hash(167726968), k, seqId, 0, false),
    //     PacBio::Pancake::SeedDB::Seed::Encode(Hash(189282306), k, seqId, 2, true),
    //     PacBio::Pancake::SeedDB::Seed::Encode(Hash(315756032), k, seqId, 3, true),
    //     PacBio::Pancake::SeedDB::Seed::Encode(Hash(259553542), k, seqId, 8, false),
    //     PacBio::Pancake::SeedDB::Seed::Encode(Hash(364792648), k, seqId, 10, true),
    //     PacBio::Pancake::SeedDB::Seed::Encode(Hash(126904899), k, seqId, 14, false),
    //     PacBio::Pancake::SeedDB::Seed::Encode(Hash(80713150), k, seqId, 18, true),
    //     PacBio::Pancake::SeedDB::Seed::Encode(Hash(27856109), k, seqId, 19, false),
    //     PacBio::Pancake::SeedDB::Seed::Encode(Hash(76041465), k, seqId, 24, true),
    //     PacBio::Pancake::SeedDB::Seed::Encode(Hash(54228923), k, seqId, 26, false),
    //     PacBio::Pancake::SeedDB::Seed::Encode(Hash(55644705), k, seqId, 31, true),
    //     PacBio::Pancake::SeedDB::Seed::Encode(Hash(502196464), k, seqId, 33, false),
    //     PacBio::Pancake::SeedDB::Seed::Encode(Hash(518950656), k, seqId, 35, false),
    //     PacBio::Pancake::SeedDB::Seed::Encode(Hash(536855377), k, seqId, 39, true),
    //     PacBio::Pancake::SeedDB::Seed::Encode(Hash(503315509), k, seqId, 41, true),
    //     PacBio::Pancake::SeedDB::Seed::Encode(Hash(125828877), k, seqId, 42, true),
    //     PacBio::Pancake::SeedDB::Seed::Encode(Hash(74973168), k, seqId, 44, true),
    //     PacBio::Pancake::SeedDB::Seed::Encode(Hash(35767), k, seqId, 46, false),
    //     PacBio::Pancake::SeedDB::Seed::Encode(Hash(143070), k, seqId, 47, false),
    //     PacBio::Pancake::SeedDB::Seed::Encode(Hash(572280), k, seqId, 48, false),
    //     PacBio::Pancake::SeedDB::Seed::Encode(Hash(2289123), k, seqId, 49, false),
    //     PacBio::Pancake::SeedDB::Seed::Encode(Hash(9156492), k, seqId, 50, false),
    //     PacBio::Pancake::SeedDB::Seed::Encode(Hash(36625970), k, seqId, 51, false),
    // };

    /*
     * WITH the invertible hash function, these would be the generated seeds.
    */
    const std::vector<PacBio::Pancake::Int128t> expectedSeeds = {
        PacBio::Pancake::SeedDB::Seed::Encode(132144311, 15, 123, 4, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(33567509, 15, 123, 5, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(18410160, 15, 123, 6, 0),
        PacBio::Pancake::SeedDB::Seed::Encode(457467302, 15, 123, 10, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(60328715, 15, 123, 12, 0),
        PacBio::Pancake::SeedDB::Seed::Encode(35576091, 15, 123, 15, 0),
        PacBio::Pancake::SeedDB::Seed::Encode(65050698, 15, 123, 20, 0),
        PacBio::Pancake::SeedDB::Seed::Encode(107140825, 15, 123, 23, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(330443040, 15, 123, 24, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(668958223, 15, 123, 27, 0),
        PacBio::Pancake::SeedDB::Seed::Encode(152777239, 15, 123, 30, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(235365058, 15, 123, 34, 0),
        PacBio::Pancake::SeedDB::Seed::Encode(39491807, 15, 123, 37, 0),
        PacBio::Pancake::SeedDB::Seed::Encode(156927623, 15, 123, 38, 0),
        PacBio::Pancake::SeedDB::Seed::Encode(511189515, 15, 123, 41, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(394814865, 15, 123, 46, 0),
        PacBio::Pancake::SeedDB::Seed::Encode(222756012, 15, 123, 47, 0),
        PacBio::Pancake::SeedDB::Seed::Encode(423331484, 15, 123, 52, 0),
        PacBio::Pancake::SeedDB::Seed::Encode(70557206, 15, 123, 54, 0),
    };

    HelperTestGenerateMinimizers(seq, seqId, k, w, space, useHPC, maxHPCLen, useRC, expectedRv,
                                 expectedSeeds, expectedThrow);
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
    const bool expectedThrow = false;
    /*
     * Without the invertible hash function, these would be the generated seeds.
    */
    // const std::vector<PacBio::Pancake::Int128t> expectedSeeds = {
    //     PacBio::Pancake::SeedDB::Seed::Encode(Hash(2), k, seqId, 2, true),
    //     PacBio::Pancake::SeedDB::Seed::Encode(Hash(56), k, seqId, 6, true),
    //     PacBio::Pancake::SeedDB::Seed::Encode(Hash(131), k, seqId, 8, true),
    //     PacBio::Pancake::SeedDB::Seed::Encode(Hash(180), k, seqId, 12, true),
    //     PacBio::Pancake::SeedDB::Seed::Encode(Hash(121), k, seqId, 14, false),
    // };
    /*
     * WITH the invertible hash function, these would be the generated seeds.
    */
    const std::vector<PacBio::Pancake::Int128t> expectedSeeds = {
        PacBio::Pancake::SeedDB::Seed::Encode(67, 5, 123, 3, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(195, 5, 123, 4, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(227, 5, 123, 5, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(235, 5, 123, 6, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(419, 5, 123, 9, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(351, 5, 123, 12, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(239, 5, 123, 15, 0),
    };

    HelperTestGenerateMinimizers(seq, seqId, k, w, space, useHPC, maxHPCLen, useRC, expectedRv,
                                 expectedSeeds, expectedThrow);
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
    const bool expectedThrow = false;
    /*
     * Without the invertible hash function, these would be the generated seeds.
    */
    // const std::vector<PacBio::Pancake::Int128t> expectedSeeds = {
    //     PacBio::Pancake::SeedDB::Seed::Encode(2, seqId, 2, true),
    //     PacBio::Pancake::SeedDB::Seed::Encode(56, seqId, 6, true),
    //     PacBio::Pancake::SeedDB::Seed::Encode(131, seqId, 8, true),
    //     PacBio::Pancake::SeedDB::Seed::Encode(180, seqId, 12, true),
    //     PacBio::Pancake::SeedDB::Seed::Encode(121, seqId, 14, false),
    //     PacBio::Pancake::SeedDB::Seed::Encode(2, seqId, 33, true),
    //     PacBio::Pancake::SeedDB::Seed::Encode(56, seqId, 37, true),
    //     PacBio::Pancake::SeedDB::Seed::Encode(131, seqId, 39, true),
    //     PacBio::Pancake::SeedDB::Seed::Encode(180, seqId, 43, true),
    //     PacBio::Pancake::SeedDB::Seed::Encode(121, seqId, 45, false),
    // };
    /*
     * WITH the invertible hash function, these would be the generated seeds.
    */
    const std::vector<PacBio::Pancake::Int128t> expectedSeeds = {
        PacBio::Pancake::SeedDB::Seed::Encode(67, 5, 123, 3, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(195, 5, 123, 4, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(227, 5, 123, 5, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(235, 5, 123, 6, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(419, 5, 123, 9, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(351, 5, 123, 12, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(239, 5, 123, 15, 0),
        PacBio::Pancake::SeedDB::Seed::Encode(67, 5, 123, 34, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(195, 5, 123, 35, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(227, 5, 123, 36, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(235, 5, 123, 37, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(419, 5, 123, 40, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(351, 5, 123, 43, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(239, 5, 123, 46, 0),
    };

    HelperTestGenerateMinimizers(seq, seqId, k, w, space, useHPC, maxHPCLen, useRC, expectedRv,
                                 expectedSeeds, expectedThrow);
}

TEST(GenerateMinimizers, SeedSize32BasePairs_PolyA)
{
    // Inputs.
    const std::string seq = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
    const int32_t seqId = 0;
    const int32_t k = 32;
    const int32_t w = 1;
    const int32_t space = 0;
    const bool useHPC = false;
    const int32_t maxHPCLen = 10;
    const bool useRC = true;

    // Expected results.
    const int32_t expectedRv = 0;
    const bool expectedThrow = true;
    const std::vector<PacBio::Pancake::Int128t> expectedSeeds = {};

    HelperTestGenerateMinimizers(seq, seqId, k, w, space, useHPC, maxHPCLen, useRC, expectedRv,
                                 expectedSeeds, expectedThrow);
}

TEST(GenerateMinimizers, SeedSize28BasePairs_PolyA)
{
    // Inputs.
    const std::string seq = "AAAAAAAAAAAAAAAAAAAAAAAAAAAA";
    const int32_t seqId = 0;
    const int32_t k = 28;
    const int32_t w = 1;
    const int32_t space = 0;
    const bool useHPC = false;
    const int32_t maxHPCLen = 10;
    const bool useRC = true;

    // Helper function.
    const uint64_t mask = SeedDB::ComputeKmerMask(k);
    auto Hash = [&](uint64_t val) { return SeedDB::InvertibleHash(val, mask); };

    // Expected results.
    const int32_t expectedRv = 0;
    const bool expectedThrow = false;
    const std::vector<PacBio::Pancake::Int128t> expectedSeeds = {
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(0), k, seqId, 0, false),
    };

    HelperTestGenerateMinimizers(seq, seqId, k, w, space, useHPC, maxHPCLen, useRC, expectedRv,
                                 expectedSeeds, expectedThrow);
}

TEST(GenerateMinimizers, SeedSize28BasePairs_PolyT_WithRC)
{
    // Inputs.
    const std::string seq = "TTTTTTTTTTTTTTTTTTTTTTTTTTTT";
    const int32_t seqId = 0;
    const int32_t k = 28;
    const int32_t w = 1;
    const int32_t space = 0;
    const bool useHPC = false;
    const int32_t maxHPCLen = 10;
    const bool useRC = true;

    // Helper function.
    const uint64_t mask = SeedDB::ComputeKmerMask(k);
    auto Hash = [&](uint64_t val) { return SeedDB::InvertibleHash(val, mask); };

    // Expected results.
    const int32_t expectedRv = 0;
    const bool expectedThrow = false;
    const std::vector<PacBio::Pancake::Int128t> expectedSeeds = {
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(0), k, seqId, 0, true),
    };

    HelperTestGenerateMinimizers(seq, seqId, k, w, space, useHPC, maxHPCLen, useRC, expectedRv,
                                 expectedSeeds, expectedThrow);
}

TEST(GenerateMinimizers, SeedSize28BasePairs_PolyT_OnlyFWD)
{
    // Inputs.
    const std::string seq = "TTTTTTTTTTTTTTTTTTTTTTTTTTTT";
    const int32_t seqId = 0;
    const int32_t k = 28;
    const int32_t w = 1;
    const int32_t space = 0;
    const bool useHPC = false;
    const int32_t maxHPCLen = 10;
    const bool useRC = false;

    // Helper function.
    const uint64_t mask = SeedDB::ComputeKmerMask(k);
    auto Hash = [&](uint64_t val) { return SeedDB::InvertibleHash(val, mask); };

    // Expected results.
    const int32_t expectedRv = 0;
    const bool expectedThrow = false;
    const std::vector<PacBio::Pancake::Int128t> expectedSeeds = {
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(0xFFFFFFFFFFFFFFFF), k, seqId, 0, false),
    };

    HelperTestGenerateMinimizers(seq, seqId, k, w, space, useHPC, maxHPCLen, useRC, expectedRv,
                                 expectedSeeds, expectedThrow);
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

    // Helper function.
    const uint64_t mask = SeedDB::ComputeKmerMask(k);
    auto Hash = [&](uint64_t val) { return SeedDB::InvertibleHash(val, mask); };

    // Expected results.
    const int32_t expectedRv = 0;
    const bool expectedThrow = false;
    const std::vector<PacBio::Pancake::Int128t> expectedSeeds = {
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(110), 8, seqId, 0, false),
    };

    HelperTestGenerateMinimizers(seq, seqId, k, w, space, useHPC, maxHPCLen, useRC, expectedRv,
                                 expectedSeeds, expectedThrow);
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
    std::vector<PacBio::Pancake::Int128t> seedsWithHP;
    int32_t rvWithHP = 0;
    {
        const uint8_t* seqData = reinterpret_cast<const uint8_t*>(seqWithHP.data());
        int32_t seqLen = seqWithHP.size();
        rvWithHP = PacBio::Pancake::SeedDB::GenerateMinimizers(
            seedsWithHP, seqData, seqLen, 0, seqId, k, w, space, useRC, true, maxHPCLen);

        // Reset the span because this will differ between the no-hpc and with-hpc runs.
        for (auto& val : seedsWithHP) {
            auto s = PacBio::Pancake::SeedDB::Seed(val);
            s.span = 0;
            val = s.To128t();
        }
    }

    // Run the unit under test WITHOUT homopolymer compression on the sequence with no homopolymers.
    std::vector<PacBio::Pancake::Int128t> seedsNoHP;
    int32_t rvNoHP = 0;
    {
        const uint8_t* seqData = reinterpret_cast<const uint8_t*>(seqWithoutHP.data());
        int32_t seqLen = seqWithoutHP.size();
        rvNoHP = PacBio::Pancake::SeedDB::GenerateMinimizers(seedsNoHP, seqData, seqLen, 0, seqId,
                                                             k, w, space, useRC, false, maxHPCLen);

        // Reset the span because this will differ between the no-hpc and with-hpc runs.
        for (auto& val : seedsNoHP) {
            auto s = PacBio::Pancake::SeedDB::Seed(val);
            s.span = 0;
            val = s.To128t();
        }
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
    }
}

TEST(GenerateMinimizers, SpacedSeed_Space1_31bp_JustOneSeed)
{
    /*
     * Here we test the spaced seed construction, with 1 skipped base
     * in between every 2 inclusive bases (i.e. space = 1).
     * This should skip every 'T' base and leave only 'A' bases in the test.
    */
    // Inputs.
    const std::string seq = "ATATATATATATATATATATATATATATATA";
    const int32_t seqId = 0;
    const int32_t k = 16;
    const int32_t w = 1;
    const int32_t space = 1;
    const bool useHPC = false;
    const int32_t maxHPCLen = 10;
    const bool useRC = false;

    // Helper function.
    const uint64_t mask = SeedDB::ComputeKmerMask(k);
    auto Hash = [&](uint64_t val) { return SeedDB::InvertibleHash(val, mask); };

    // Expected results.
    const int32_t expectedRv = 0;
    const bool expectedThrow = false;
    const std::vector<PacBio::Pancake::Int128t> expectedSeeds = {
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(0), 31, seqId, 0, false),
    };

    HelperTestGenerateMinimizers(seq, seqId, k, w, space, useHPC, maxHPCLen, useRC, expectedRv,
                                 expectedSeeds, expectedThrow);
}

TEST(GenerateMinimizers, SpacedSeed_Space1_32bp_TwoSeeds)
{
    /*
     * Here we test the spaced seed construction, with 1 skipped base
     * in between every 2 inclusive bases (i.e. space = 1).
     * This should skip every 'T' base and leave only 'A' bases in the test.
    */
    // Inputs.
    const std::string seq = "TATATATATATATATATATATATATATATATA";
    const int32_t seqId = 0;
    const int32_t k = 16;
    const int32_t w = 1;
    const int32_t space = 1;
    const bool useHPC = false;
    const int32_t maxHPCLen = 10;
    const bool useRC = false;

    // Helper function.
    const uint64_t mask = SeedDB::ComputeKmerMask(k);
    auto Hash = [&](uint64_t val) { return SeedDB::InvertibleHash(val, mask); };

    // Expected results.
    const int32_t expectedRv = 0;
    const bool expectedThrow = false;
    const std::vector<PacBio::Pancake::Int128t> expectedSeeds = {
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(0x0FFFFFFFF), 31, seqId, 0,
                                              false),  // 16 bases of 'T's.
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(0), 31, seqId, 1, false),
    };

    HelperTestGenerateMinimizers(seq, seqId, k, w, space, useHPC, maxHPCLen, useRC, expectedRv,
                                 expectedSeeds, expectedThrow);
}

TEST(GenerateMinimizers, SpacedSeed_Space1_33bp_ThreeSeeds)
{
    /*
     * Here we test the spaced seed construction, with 1 skipped base
     * in between every 2 inclusive bases (i.e. space = 1).
     * This should skip every 'T' base and leave only 'A' bases in the test.
    */
    // Inputs.
    const std::string seq = "TATATATATATATATATATATATATATATATAG";
    const int32_t seqId = 0;
    const int32_t k = 16;
    const int32_t w = 1;
    const int32_t space = 1;
    const bool useHPC = false;
    const int32_t maxHPCLen = 10;
    const bool useRC = false;

    // Helper function.
    const uint64_t mask = SeedDB::ComputeKmerMask(k);
    auto Hash = [&](uint64_t val) { return SeedDB::InvertibleHash(val, mask); };

    // Expected results.
    const int32_t expectedRv = 0;
    const bool expectedThrow = false;
    const std::vector<PacBio::Pancake::Int128t> expectedSeeds = {
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(0x0FFFFFFFF), 31, seqId, 0,
                                              false),  // 16 bases of 'T's.
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(0), 31, seqId, 1, false),
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(0x0FFFFFFFE), 31, seqId, 2, false),
    };

    HelperTestGenerateMinimizers(seq, seqId, k, w, space, useHPC, maxHPCLen, useRC, expectedRv,
                                 expectedSeeds, expectedThrow);
}

TEST(GenerateMinimizers, FromStrings)
{
    // Inputs.
    const int32_t k = 5;
    const int32_t w = 4;
    const int32_t space = 0;
    const bool useHPC = false;
    const int32_t maxHPCLen = 10;
    const bool useRC = true;
    const std::vector<std::string> seqs = {
        "AAAAAAAAAA", "AGCTTTTCATTCTGACTGCANNNACTNNNNNAGCTTTTCATTCTGACTGCA",
    };

    // Helper function.
    const uint64_t mask = SeedDB::ComputeKmerMask(k);
    auto Hash = [&](uint64_t val) { return SeedDB::InvertibleHash(val, mask); };

    // Expected results.
    const std::vector<PacBio::Pancake::Int128t> expectedSeeds = {
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(0), k, 0, 3, false),
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(0), k, 0, 4, false),
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(0), k, 0, 5, false),

        PacBio::Pancake::SeedDB::Seed::Encode(67, 5, 1, 3, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(195, 5, 1, 4, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(227, 5, 1, 5, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(235, 5, 1, 6, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(419, 5, 1, 9, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(351, 5, 1, 12, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(239, 5, 1, 15, 0),
        PacBio::Pancake::SeedDB::Seed::Encode(67, 5, 1, 34, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(195, 5, 1, 35, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(227, 5, 1, 36, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(235, 5, 1, 37, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(419, 5, 1, 40, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(351, 5, 1, 43, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(239, 5, 1, 46, 0),
    };
    std::vector<int32_t> expectedSequenceLengths = {10, 51};

    // Run unit under test.
    std::vector<PacBio::Pancake::Int128t> results;
    std::vector<int32_t> sequenceLengths;
    PacBio::Pancake::SeedDB::GenerateMinimizers(results, sequenceLengths, seqs, k, w, space, useRC,
                                                useHPC, maxHPCLen);

    // Evaluate.
    EXPECT_EQ(expectedSeeds, results);
    EXPECT_EQ(expectedSequenceLengths, sequenceLengths);
}
