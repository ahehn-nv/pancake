#include <gtest/gtest.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <tuple>
#include <vector>

#include <pacbio/alignment/AlignmentTools.h>
#include <pbcopper/third-party/edlib.h>

#include "PancakeTestData.h"

using namespace PacBio::Pancake;

namespace AlignmentToolsTests {

TEST(Test_AlignmentTools, EmptyInput)
{
    std::vector<unsigned char> input = {};
    PacBio::BAM::Cigar expected;
    PacBio::BAM::Cigar result = EdlibAlignmentToCigar(input.data(), input.size());
    EXPECT_EQ(expected, result);
}

TEST(Test_AlignmentTools, SimpleTest)
{
    // Edlib move codes: 0: '=', 1: 'I', 2: 'D', 3: 'X'
    std::vector<unsigned char> input = {EDLIB_EDOP_MATCH,    EDLIB_EDOP_MATCH,  EDLIB_EDOP_MATCH,
                                        EDLIB_EDOP_MISMATCH, EDLIB_EDOP_INSERT, EDLIB_EDOP_DELETE,
                                        EDLIB_EDOP_DELETE,   EDLIB_EDOP_INSERT};
    PacBio::BAM::Cigar expected("3=1X1I2D1I");
    PacBio::BAM::Cigar result = EdlibAlignmentToCigar(input.data(), input.size());
    EXPECT_EQ(expected, result);
}

TEST(Test_AlignmentTools_ValidateCigar, ArrayOfTests)
{
    // clang-format off
    std::vector<std::tuple<std::string, std::string, std::string, std::string, bool>> testData = {
        // Good cases.
        {"Empty input", "", "", "", false},
        {"Exact match, good", "ACTG", "ACTG", "4=", false},
        {"Cigar with indels", "ACTT", "ACTG", "3=1I1D", false},
        {"Soft clipping", "AAAACTTAAA", "ACTG", "3S3=1X3S", false},
        {"Hard clipping", "ACTT", "ACTG", "3H3=1X3H", false},
        {"Simple test case", "GGATCAGTTTTATATACAC", "GAGTTCGTTCTATATATACAC", "1=1I1=1D1=1D1=1I3=1X3=2D4=2I", false},

        // Bad cases.
        {"One base is incorrect", "ACTT", "ACTG", "4=", true},
        {"Length mismatch", "ACT", "ACTG", "4=", true},
        {"Non-zero CIGAR, but empty seqs", "", "", "4=", true},
        {"Length mismatch with indels", "ACTT", "ACTG", "3=1I1D1I", true},
        {"Length mismatch with mismatches", "ACTT", "ACTG", "3=2X", true},

    };
    // clang-format on

    for (const auto& data : testData) {
        // Inputs.
        const std::string testName = std::get<0>(data);
        const std::string& query = std::get<1>(data);
        const std::string& target = std::get<2>(data);
        PacBio::BAM::Cigar cigar(std::get<3>(data));
        bool shouldThrow = std::get<4>(data);

        // Name the test.
        SCOPED_TRACE("ValidateCigar-" + testName);

        // Run.
        if (shouldThrow) {
            EXPECT_THROW(
                {
                    PacBio::Pancake::ValidateCigar(query.c_str(), query.size(), target.c_str(),
                                                   target.size(), cigar, "");
                },
                std::runtime_error);
        } else {
            PacBio::Pancake::ValidateCigar(query.c_str(), query.size(), target.c_str(),
                                           target.size(), cigar, "");
        }
    }
}

class ExtractVariantStringTestCase
{
public:
    std::string testName;
    std::string query;
    std::string target;
    std::string cigar;
    bool maskHomopolymers;
    bool maskSimpleRepeats;
    std::string expectedQueryVariants;
    std::string expectedTargetVariants;
    PacBio::Pancake::Alignment::DiffCounts expectedDiffsPerBase;
    PacBio::Pancake::Alignment::DiffCounts expectedDiffsPerEvent;
};

TEST(Test_AlignmentTools_ExtractVariantString, ArrayOfTests)
{
    // clang-format off
    std::vector<ExtractVariantStringTestCase> testData = {
        ExtractVariantStringTestCase{"Empty input", "", "", "", false, false, "", "", {}, {}},

        /*
        * Q: GGA-T-CAGTTTT AT--ATACAC
        *    | | | | |||X| ||  ||||
        * T: G-AGTTC-GTTCT ATATATAC--
        * Query: indel G is in a homopolymer, A is not in a HP.
        * Target: indel G is not a HP, and T is.
        * There is one insertion event of 2bp and one deletion event of 2bp. Both are in simple repeats.
        */
        ExtractVariantStringTestCase{"Simple test case",
                                     "GGATCAGTTTTATATACAC", "GAGTTCGTTCTATATATACAC", "1=1I1=1D1=1D1=1I3=1X3=2D4=2I", false, false, "GATAC", "GTCAT", {14, 1, 4, 4}, {14, 1, 3, 3}},

        ExtractVariantStringTestCase{"Test case with HP masking",
                                     "GGATCAGTTTTATATACAC", "GAGTTCGTTCTATATATACAC", "1=1I1=1D1=1D1=1I3=1X3=2D4=2I", true, false, "gATAC", "GtCAT", {14, 1, 3, 3}, {14, 1, 2, 2}},

        ExtractVariantStringTestCase{"Test case with simple repeat masking",
                                     "GGATCAGTTTTATATACAC", "GAGTTCGTTCTATATATACAC", "1=1I1=1D1=1D1=1I3=1X3=2D4=2I", false, true, "GATac", "GTCat", {14, 1, 2, 2}, {14, 1, 2, 2}},

        ExtractVariantStringTestCase{"Test case with HP and simple repeat masking",
                                     "GGATCAGTTTTATATACAC", "GAGTTCGTTCTATATATACAC", "1=1I1=1D1=1D1=1I3=1X3=2D4=2I", true, true, "gATac", "GtCat", {14, 1, 1, 1}, {14, 1, 1, 1}},

        /*
         * This test case ("GCAC', "GAC") detected a false masking bug in the code. The HP masking used
         * to check forward in target for insertions and backward in target for deletions. This is wrong
         * because the coordinate in the strand with the gap points to the base _after_ the gap (and not before
         * like I originally thought). So it should be vice versa.
        */
        ExtractVariantStringTestCase{"Small HP masking test 1 - insertion, no masking",
                                     "GCAC", "GAC", "1=1I2=", true, true, "C", "", {3, 0, 1, 0}, {3, 0, 1, 0}},
        ExtractVariantStringTestCase{"Small HP masking test 1 - deletion, no masking",
                                     "GAC", "GCAC", "1=1D2=", true, true, "", "C", {3, 0, 0, 1}, {3, 0, 0, 1}},

        ExtractVariantStringTestCase{"Small HP masking test 2 - insertion, no masking",
                                     "TGCC", "TCC", "1=1I2=", true, true, "G", "", {3, 0, 1, 0}, {3, 0, 1, 0}},
        ExtractVariantStringTestCase{"Small HP masking test 2 - deletion, no masking",
                                     "TCC", "TGCC", "1=1D2=", true, true, "", "G", {3, 0, 0, 1}, {3, 0, 0, 1}},

        ExtractVariantStringTestCase{"Small HP masking test 3 - insertion, no masking",
                                     "ACAA", "AAA", "1=1I2=", true, true, "C", "", {3, 0, 1, 0}, {3, 0, 1, 0}},
        ExtractVariantStringTestCase{"Small HP masking test 3 - insertion, no masking",
                                     "AAA", "ACAA", "1=1D2=", true, true, "", "C", {3, 0, 0, 1}, {3, 0, 0, 1}},


    };
    // clang-format on

    for (const auto& data : testData) {
        // Name the test.
        SCOPED_TRACE("ExtractVariantString-" + data.testName);
        PacBio::BAM::Cigar cigar(data.cigar);
        std::string resultQueryVariants;
        std::string resultTargetVariants;
        PacBio::Pancake::Alignment::DiffCounts resultDiffsPerBase;
        PacBio::Pancake::Alignment::DiffCounts resultDiffsPerEvent;

        // Run.
        PacBio::Pancake::ExtractVariantString(
            data.query.c_str(), data.query.size(), data.target.c_str(), data.target.size(), cigar,
            data.maskHomopolymers, data.maskSimpleRepeats, resultQueryVariants,
            resultTargetVariants, resultDiffsPerBase, resultDiffsPerEvent);

        // Evaluate.
        EXPECT_EQ(data.expectedQueryVariants, resultQueryVariants);
        EXPECT_EQ(data.expectedTargetVariants, resultTargetVariants);
        EXPECT_EQ(data.expectedDiffsPerBase, resultDiffsPerBase);
        EXPECT_EQ(data.expectedDiffsPerEvent, resultDiffsPerEvent);
    }
}

TEST(Test_AlignmentTools_FindTargetPosFromCigar, ArrayOfTests)
{
    // clang-format off
    std::vector<std::tuple<std::string, std::string, int32_t, int32_t, bool>> testData = {
        {"Empty input", "", 0, -1, true},
        {"Simple matches", "10=", 4, 4, false},
        {"Simple matches", "10=", -1, -1, true},
        {"Simple matches", "10=", 10, -1, true},

        {"Simple CIGAR pos 0", "2=2D1=1X1=2I1=", 0, 0, false},
        {"Simple CIGAR pos 1", "2=2D1=1X1=2I1=", 1, 1, false},
        {"Simple CIGAR pos 2", "2=2D1=1X1=2I1=", 2, 4, false},
        {"Simple CIGAR pos 3", "2=2D1=1X1=2I1=", 3, 5, false},
        {"Simple CIGAR pos 4", "2=2D1=1X1=2I1=", 4, 6, false},
        {"Simple CIGAR pos 5", "2=2D1=1X1=2I1=", 5, 6, false},
        {"Simple CIGAR pos 6", "2=2D1=1X1=2I1=", 6, 6, false},
        {"Simple CIGAR pos 7", "2=2D1=1X1=2I1=", 7, 7, false},

        {"Custom real 1", "59=1I10=1I601=1I48=1D274=1D432=1I84=1I573=1I94=1D1545=1I1=1D1131=1D1042=1I329=1I581=1D452=", 1380, 1379, false},
    };
    // clang-format on

    for (const auto& data : testData) {
        // Get the data.
        const std::string testName = std::get<0>(data);
        PacBio::BAM::Cigar cigar(std::get<1>(data));
        int32_t queryPos = std::get<2>(data);
        int32_t expected = std::get<3>(data);
        bool shouldThrow = std::get<4>(data);

        // Name the test.
        SCOPED_TRACE("FindTargetPosFromCigar-" + testName);

        // Run.
        if (shouldThrow) {
            EXPECT_THROW({ PacBio::Pancake::FindTargetPosFromCigar(cigar, queryPos); },
                         std::runtime_error);
        } else {
            int32_t result = PacBio::Pancake::FindTargetPosFromCigar(cigar, queryPos);
            // Evaluate.
            EXPECT_EQ(expected, result);
        }
    }
}

TEST(Test_AlignmentTools_ComputeDiffCounts, ArrayOfTests)
{
    // clang-format off
    std::vector<std::tuple<std::string, std::string, std::string, std::string, Alignment::DiffCounts, bool>> testData = {
        {"Empty input", "", "", "", Alignment::DiffCounts(), false},
        {"Exact match", "4=", "", "", Alignment::DiffCounts(4, 0, 0, 0), false},
        {"Variant pos out of bounds (X)", "3=1X1=", "", "", Alignment::DiffCounts(), true},
        {"Variant pos out of bounds (I)", "3=1I1=", "", "", Alignment::DiffCounts(), true},
        {"Variant pos out of bounds (D)", "3=1D1=", "", "", Alignment::DiffCounts(), true},
        {"Valid case, 1X", "3=1X1=", "A", "C", Alignment::DiffCounts(4, 1, 0, 0), false},
        {"Valid case, 1X, masked", "3=1X1=", "a", "c", Alignment::DiffCounts(4, 0, 0, 0), false},
        {"Valid case, 1X, mixed mask, should throw", "3=1X1=", "a", "C", Alignment::DiffCounts(), true},
        {"Valid case with multiple diffs and 1 masked I", "3=1X1=1D2=2I2=1I", "ATAt", "CT", Alignment::DiffCounts(8, 1, 2, 1), false},
        {"Valid case with clipping and reference skip", "3S3=1X1=1D2=1I2N2=1I3H", "ATt", "CT", Alignment::DiffCounts(8, 1, 1, 1), false},
        {"Bad case, the 2-base I has mixed masking", "3=1X1=1D2=2I2=1I", "ATat", "CT", Alignment::DiffCounts(), true},
        {"Bad case, the 2-base B has mixed masking", "3=1X1=2D2=2I2=1I", "ATat", "CTc", Alignment::DiffCounts(), true},
    };
    // clang-format on

    for (const auto& data : testData) {
        // Inputs.
        const std::string testName = std::get<0>(data);
        PacBio::BAM::Cigar cigar(std::get<1>(data));
        const std::string& queryVariants = std::get<2>(data);
        const std::string& targetVariants = std::get<3>(data);
        const PacBio::Pancake::Alignment::DiffCounts& expected = std::get<4>(data);
        bool shouldThrow = std::get<5>(data);

        // Name the test.
        SCOPED_TRACE("ComputeDiffCounts-" + testName);

        // Run.
        if (shouldThrow) {
            EXPECT_THROW(
                { PacBio::Pancake::ComputeDiffCounts(cigar, queryVariants, targetVariants); },
                std::runtime_error);
        } else {
            auto result = PacBio::Pancake::ComputeDiffCounts(cigar, queryVariants, targetVariants);
            // Evaluate.
            EXPECT_EQ(expected, result);
        }
    }
}
}
