#include <gtest/gtest.h>
#include <pacbio/seqdb/Twobit.h>

void HelperRoundTrip(const std::string& inBases, int32_t numBases, const std::vector<uint8_t>& expectedTwobit,
                     const std::vector<PacBio::Pancake::Range>& expectedRanges,
                     bool expectedException)
{
    // Compress the input sequence.
    std::vector<uint8_t> twobit;
    std::vector<PacBio::Pancake::Range> ranges;
    PacBio::Pancake::CompressSequence(inBases, twobit, ranges);

    // Check the compression results.
    EXPECT_EQ(expectedTwobit, twobit);
    EXPECT_EQ(expectedRanges, ranges);

    // Decompress the sequence.
    // Run round trip, expected to see the same sequence or an exception.
    bool exceptionCaught = false;
    try {
        std::string roundResult;
        PacBio::Pancake::DecompressSequence(twobit, numBases, ranges, roundResult);
        EXPECT_EQ(inBases, roundResult);
    } catch (const std::runtime_error& e) {
        exceptionCaught = true;
    }
    EXPECT_EQ(expectedException, exceptionCaught);
}

TEST(Twobit_RoundTrip, Empty)
{
    const std::string inBases("");
    const int32_t numBases = inBases.size();
    const std::vector<uint8_t> expectedTwobit;
    const std::vector<PacBio::Pancake::Range> expectedRanges;
    HelperRoundTrip(inBases, numBases, expectedTwobit, expectedRanges, false);
}

TEST(Twobit_RoundTrip, Simple_1_Base_NoN_1)
{
    // 'A' = 00-00-00-00(binary) = 0(dec)
    const std::string inBases("A");
    const int32_t numBases = inBases.size();
    const std::vector<uint8_t> expectedTwobit{0};
    const std::vector<PacBio::Pancake::Range> expectedRanges{{0, 1}};
    HelperRoundTrip(inBases, numBases, expectedTwobit, expectedRanges, false);
}

TEST(Twobit_RoundTrip, Simple_1_Base_NoN_2)
{
    // 'G' = 10-00-00-00(binary) = 128(dec)
    const std::string inBases("G");
    const int32_t numBases = inBases.size();
    const std::vector<uint8_t> expectedTwobit{128};
    const std::vector<PacBio::Pancake::Range> expectedRanges{{0, 1}};
    HelperRoundTrip(inBases, numBases, expectedTwobit, expectedRanges, false);
}

TEST(Twobit_RoundTrip, Simple_2_Bases_NoN)
{
    // 'GC' = 10-01-00-00(binary) = 144(dec)
    const std::string inBases("GC");
    const int32_t numBases = inBases.size();
    const std::vector<uint8_t> expectedTwobit{144};
    const std::vector<PacBio::Pancake::Range> expectedRanges{{0, 2}};
    HelperRoundTrip(inBases, numBases, expectedTwobit, expectedRanges, false);
}

TEST(Twobit_RoundTrip, Simple_3_Bases_NoN)
{
    // 'GCT' = 10-01-11-00(binary) = 156(dec)
    const std::string inBases("GCT");
    const int32_t numBases = inBases.size();
    const std::vector<uint8_t> expectedTwobit{156};
    const std::vector<PacBio::Pancake::Range> expectedRanges{{0, 3}};
    HelperRoundTrip(inBases, numBases, expectedTwobit, expectedRanges, false);
}

TEST(Twobit_RoundTrip, Simple_4_Bases_NoN)
{
    // ACGT = 00-01-10-11(binary) = 27(dec)
    const std::string inBases("ACGT");
    const int32_t numBases = inBases.size();
    const std::vector<uint8_t> expectedTwobit{27};
    const std::vector<PacBio::Pancake::Range> expectedRanges{{0, 4}};
    HelperRoundTrip(inBases, numBases, expectedTwobit, expectedRanges, false);
}

TEST(Twobit_RoundTrip, Simple_WithN_1)
{
    // ACGTNTTT = 00-01-10-11(binary) 11-11-11-00(binary) = 27(dec) 252(dec)
    const std::string inBases("ACGTNTTT");
    const int32_t numBases = inBases.size();
    const std::vector<uint8_t> expectedTwobit{27, 252};
    const std::vector<PacBio::Pancake::Range> expectedRanges{{0, 4}, {5, 8}};
    HelperRoundTrip(inBases, numBases, expectedTwobit, expectedRanges, false);
}

TEST(Twobit_RoundTrip, Simple_WithN_2)
{
    // GCTT AT = 10-01-11-11(binary) 00-11-00-00(binary) = 159(dec) 48(dec)
    const std::string inBases("GCTNTAT");
    const int32_t numBases = inBases.size();
    const std::vector<uint8_t> expectedTwobit{159, 48};
    const std::vector<PacBio::Pancake::Range> expectedRanges{{0, 3}, {4, 7}};
    HelperRoundTrip(inBases, numBases, expectedTwobit, expectedRanges, false);
}

TEST(Twobit_RoundTrip, Simple_WithMultipleN)
{
    // GCTA CGTA T = 10-01-11-00(bin) 01-10-11-00(bin) 11-00-00-00(bin) = 156(dec) 108(dec) 192(dec)
    const std::string inBases("GCTNANCNNNGNTAT");
    const int32_t numBases = inBases.size();
    const std::vector<uint8_t> expectedTwobit{156, 108, 192};
    const std::vector<PacBio::Pancake::Range> expectedRanges{
        {0, 3}, {4, 5}, {6, 7}, {10, 11}, {12, 15}};
    HelperRoundTrip(inBases, numBases, expectedTwobit, expectedRanges, false);
}

TEST(Twobit_RoundTrip, Simple_WithN_In_Last_Pos)
{
    // GCAT = 10-01-00-11(bin) = 147(dec)
    const std::string inBases("GCATNNN");
    const int32_t numBases = inBases.size();
    const std::vector<uint8_t> expectedTwobit{147};
    const std::vector<PacBio::Pancake::Range> expectedRanges{{0, 4}};
    HelperRoundTrip(inBases, numBases, expectedTwobit, expectedRanges, false);
}

TEST(Twobit_RoundTrip, Simple_WithN_In_FirstPos)
{
    // GCAT = 10-01-00-11(bin) = 147(dec)
    const std::string inBases("NNNGCAT");
    const int32_t numBases = inBases.size();
    const std::vector<uint8_t> expectedTwobit{147};
    const std::vector<PacBio::Pancake::Range> expectedRanges{{3, 7}};
    HelperRoundTrip(inBases, numBases, expectedTwobit, expectedRanges, false);
}

TEST(Twobit_RoundTrip, DecompressSanityCheck1)
{
    // ACGT = 00-01-10-11(binary) = 27(dec)
    const std::string inBases("ACGT");
    const int32_t numBases = inBases.size();
    const std::vector<uint8_t> expectedTwobit{27};
    const std::vector<PacBio::Pancake::Range> expectedRanges{{0, 4}};
    HelperRoundTrip(inBases, numBases - 1, expectedTwobit, expectedRanges, true);
}

TEST(Twobit_RoundTrip, DecompressSanityCheck2)
{
    // ACGT = 00-01-10-11(binary) = 27(dec)
    const std::string inBases("ACGT");
    const int32_t numBases = inBases.size();
    const std::vector<uint8_t> expectedTwobit{27};
    const std::vector<PacBio::Pancake::Range> expectedRanges{{0, 4}};
    HelperRoundTrip(inBases, -1, expectedTwobit, expectedRanges, true);
}
