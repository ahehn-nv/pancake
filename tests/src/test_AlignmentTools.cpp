#include <gtest/gtest.h>

#include <fstream>
#include <iostream>
#include <sstream>

#include <pacbio/alignment/AlignmentTools.h>
#include <pbcopper/third-party/edlib.h>

#include "PancakeTestData.h"

using namespace PacBio::Pancake;

namespace EdlibAlignmentToolsTests {

TEST(Test_EdlibAlignmentTools, EmptyInput)
{
    std::vector<unsigned char> input = {};
    PacBio::BAM::Cigar expected;
    PacBio::BAM::Cigar result = EdlibAlignmentToCigar(input.data(), input.size());
    EXPECT_EQ(expected, result);
}

TEST(Test_EdlibAlignmentTools, SimpleTest)
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
}
