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
}
