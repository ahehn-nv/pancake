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
        ExtractVariantStringTestCase{"Small HP masking test 1 - insertion, should not masking",
                                     "GCAC", "GAC", "1=1I2=", true, true, "C", "", {3, 0, 1, 0}, {3, 0, 1, 0}},
        ExtractVariantStringTestCase{"Small HP masking test 1 - deletion, should not masking",
                                     "GAC", "GCAC", "1=1D2=", true, true, "", "C", {3, 0, 0, 1}, {3, 0, 0, 1}},
        /*
         * This is extracted around the same variant (from the dataset), but in the other direction.
         * There were no issues in this direction.
        */
        ExtractVariantStringTestCase{"Small HP masking test 2 - insertion, should not masking",
                                     "TGCC", "TCC", "1=1I2=", true, true, "G", "", {3, 0, 1, 0}, {3, 0, 1, 0}},
        ExtractVariantStringTestCase{"Small HP masking test 2 - deletion, should not masking",
                                     "TCC", "TGCC", "1=1D2=", true, true, "", "G", {3, 0, 0, 1}, {3, 0, 0, 1}},

        /*
         * Test HP insertion masking systematically.
        */
        ExtractVariantStringTestCase{"Simple HP masking test 3a - insertion, no masking",
                                     "AAACAAA", "AAAAAA", "3=1I3=", true, true, "C", "", {6, 0, 1, 0}, {6, 0, 1, 0}},
        ExtractVariantStringTestCase{"Simple HP masking test 3b - insertion, should mask",
                                     "AAACCAA", "AAACAA", "3=1I3=", true, true, "c", "", {6, 0, 0, 0}, {6, 0, 0, 0}},
        ExtractVariantStringTestCase{"Simple HP masking test 3c - insertion, no masking",
                                     "AAACACA", "AAAACA", "3=1I3=", true, true, "C", "", {6, 0, 1, 0}, {6, 0, 1, 0}},
        ExtractVariantStringTestCase{"Simple HP masking test 3d - insertion, no masking",
                                     "AAACAAC", "AAAAAC", "3=1I3=", true, true, "C", "", {6, 0, 1, 0}, {6, 0, 1, 0}},
        ExtractVariantStringTestCase{"Simple HP masking test 3e - insertion, no masking",
                                     "AACCAAA", "AACAAA", "3=1I3=", true, true, "c", "", {6, 0, 0, 0}, {6, 0, 0, 0}},
        ExtractVariantStringTestCase{"Simple HP masking test 3f - insertion, no masking",
                                     "ACACAAA", "ACAAAA", "3=1I3=", true, true, "C", "", {6, 0, 1, 0}, {6, 0, 1, 0}},
        ExtractVariantStringTestCase{"Simple HP masking test 3g - insertion, no masking",
                                     "CAACAAA", "CAAAAA", "3=1I3=", true, true, "C", "", {6, 0, 1, 0}, {6, 0, 1, 0}},
        ExtractVariantStringTestCase{"Simple HP masking test 3h - no indels, perfect match",
                                     "AAACAAA", "AAACAAA", "7=", true, true, "", "", {7, 0, 0, 0}, {7, 0, 0, 0}},
        /*
         * Test HP deletion masking systematically.
        */
        ExtractVariantStringTestCase{"Simple HP masking test 4a - deletion, no masking",
                                     "AAAAAA", "AAACAAA", "3=1D3=", true, true, "", "C", {6, 0, 0, 1}, {6, 0, 0, 1}},
        ExtractVariantStringTestCase{"Simple HP masking test 4b - deletion, should mask",
                                     "AAACAA", "AAACCAA", "3=1D3=", true, true, "", "c", {6, 0, 0, 0}, {6, 0, 0, 0}},
        ExtractVariantStringTestCase{"Simple HP masking test 4c - deletion, no masking",
                                     "AAAACA", "AAACACA", "3=1D3=", true, true, "", "C", {6, 0, 0, 1}, {6, 0, 0, 1}},
        ExtractVariantStringTestCase{"Simple HP masking test 4d - deletion, no masking",
                                     "AAAAAC", "AAACAAC", "3=1D3=", true, true, "", "C", {6, 0, 0, 1}, {6, 0, 0, 1}},
        ExtractVariantStringTestCase{"Simple HP masking test 4e - deletion, no masking",
                                     "AACAAA", "AACCAAA", "3=1D3=", true, true, "", "c", {6, 0, 0, 0}, {6, 0, 0, 0}},
        ExtractVariantStringTestCase{"Simple HP masking test 4f - deletion, no masking",
                                     "ACAAAA", "ACACAAA", "3=1D3=", true, true, "", "C", {6, 0, 0, 1}, {6, 0, 0, 1}},
        ExtractVariantStringTestCase{"Simple HP masking test 4g - deletion, no masking",
                                     "CAAAAA", "CAACAAA", "3=1D3=", true, true, "", "C", {6, 0, 0, 1}, {6, 0, 0, 1}},
        ExtractVariantStringTestCase{"Simple HP masking test 4h - no indels, perfect match",
                                     "AAACAAA", "AAACAAA", "7=", true, true, "", "", {7, 0, 0, 0}, {7, 0, 0, 0}},

        /*
         * Test simple releat insertion masking systematically.
        */
        ExtractVariantStringTestCase{"Simple repeat masking test 5a - no indel, perfect match",
                                     "AAAAAAATCAAAAAAA", "AAAAAAATCAAAAAAA", "16=", true, true, "", "", {16, 0, 0, 0}, {16, 0, 0, 0}},
        ExtractVariantStringTestCase{"Simple repeat masking test 5b - insertion, no masking",
                                     "AAAAAAATCAAAAAAA",    // AAAAAAATCAAAAAAA
                                     "AAAAAAAAAAAAAA",      // AAAAAAA--AAAAAAA
                                     "7=2I7=", true, true, "TC", "", {14, 0, 2, 0}, {14, 0, 1, 0}},
        /*
         * Sliding window for insertion masking - slide a TC event accross the target, while the query has a TC insertion.
         * Only in some cases should masking pick up.
        */
        ExtractVariantStringTestCase{"Simple repeat masking sliding window test 5a - insertion, should not mask",
                                     "AAAAAAATCAAAAAAA",    // AAAAAAATCAAAAAAA
                                     "ATCAAAAAAAAAAA",      // ATCAAAA--AAAAAAA
                                     "1=2X4=2I7=", true, true, "AATC", "TC", {12, 2, 2, 0}, {12, 2, 1, 0}},
        ExtractVariantStringTestCase{"Simple repeat masking sliding window test 5b - insertion, should not mask",
                                     "AAAAAAATCAAAAAAA",    // AAAAAAATCAAAAAAA
                                     "AATCAAAAAAAAAA",      // AATCAAA--AAAAAAA
                                     "2=2X3=2I7=", true, true, "AATC", "TC", {12, 2, 2, 0}, {12, 2, 1, 0}},
        ExtractVariantStringTestCase{"Simple repeat masking sliding window test 5c - insertion, should not mask",
                                     "AAAAAAATCAAAAAAA",    // AAAAAAATCAAAAAAA
                                     "AAATCAAAAAAAAA",      // AAATCAA--AAAAAAA
                                     "3=2X2=2I7=", true, true, "AATC", "TC", {12, 2, 2, 0}, {12, 2, 1, 0}},
        ExtractVariantStringTestCase{"Simple repeat masking sliding window test 5d - insertion, should not mask",
                                     "AAAAAAATCAAAAAAA",    // AAAAAAATCAAAAAAA
                                     "AAAATCAAAAAAAA",      // AAAATCA--AAAAAAA
                                     "4=2X1=2I7=", true, true, "AATC", "TC", {12, 2, 2, 0}, {12, 2, 1, 0}},
        ExtractVariantStringTestCase{"Simple repeat masking sliding window test 5e - insertion, should mask",
                                     "AAAAAAATCAAAAAAA",    // AAAAAAATCAAAAAAA
                                     "AAAAATCAAAAAAA",      // AAAAATC--AAAAAAA
                                     "5=2X2I7=", true, true, "AAtc", "TC", {12, 2, 0, 0}, {12, 2, 0, 0}},
        ExtractVariantStringTestCase{"Simple repeat masking sliding window test 5f - insertion, should not mask",
                                     "AAAAAAATCAAAAAAA",    // AAAAAAATCAAAAAAA
                                     "AAAAAATCAAAAAA",      // AAAAAAT--CAAAAAA
                                     "6=1X2I1X6=", true, true, "ATCA", "TC", {12, 2, 2, 0}, {12, 2, 1, 0}},
        ExtractVariantStringTestCase{"Simple repeat masking sliding window test 5g - insertion, should mask",
                                     "AAAAAAATCAAAAAAA",    // AAAAAAATCAAAAAAA
                                     "AAAAAAATCAAAAA",      // AAAAAAA--TCAAAAA
                                     "7=2I2X5=", true, true, "tcAA", "TC", {12, 2, 0, 0}, {12, 2, 0, 0}},
        ExtractVariantStringTestCase{"Simple repeat masking sliding window test 5h - insertion, should not mask",
                                     "AAAAAAATCAAAAAAA",    // AAAAAAATCAAAAAAA
                                     "AAAAAAAATCAAAA",      // AAAAAAA--ATCAAAA
                                     "7=2I1=2X4=", true, true, "TCAA", "TC", {12, 2, 2, 0}, {12, 2, 1, 0}},
        ExtractVariantStringTestCase{"Simple repeat masking sliding window test 5i - insertion, should not mask",
                                     "AAAAAAATCAAAAAAA",    // AAAAAAATCAAAAAAA
                                     "AAAAAAAAATCAAA",      // AAAAAAA--AATCAAA
                                     "7=2I2=2X3=", true, true, "TCAA", "TC", {12, 2, 2, 0}, {12, 2, 1, 0}},
        ExtractVariantStringTestCase{"Simple repeat masking sliding window test 5j - insertion, should not mask",
                                     "AAAAAAATCAAAAAAA",    // AAAAAAATCAAAAAAA
                                     "AAAAAAAAAATCAA",      // AAAAAAA--AAATCAA
                                     "7=2I3=2X2=", true, true, "TCAA", "TC", {12, 2, 2, 0}, {12, 2, 1, 0}},
        ExtractVariantStringTestCase{"Simple repeat masking sliding window test 5k - insertion, should not mask",
                                     "AAAAAAATCAAAAAAA",    // AAAAAAATCAAAAAAA
                                     "AAAAAAAAAAATCA",      // AAAAAAA--AAAATCA
                                     "7=2I4=2X1=", true, true, "TCAA", "TC", {12, 2, 2, 0}, {12, 2, 1, 0}},
        /*
         * Sliding window for deletion masking - slide a TC event accross the query, while the target has a TC insertion.
         * Only in some cases should masking pick up.
        */
        ExtractVariantStringTestCase{"Simple repeat masking sliding window test 6a - insertion, should not mask",
                                     "ATCAAAAAAAAAAA",      // ATCAAAA--AAAAAAA
                                     "AAAAAAATCAAAAAAA",    // AAAAAAATCAAAAAAA
                                     "1=2X4=2D7=", true, true, "TC", "AATC", {12, 2, 0, 2}, {12, 2, 0, 1}},
        ExtractVariantStringTestCase{"Simple repeat masking sliding window test 6b - deletion, should not mask",
                                     "AATCAAAAAAAAAA",      // AATCAAA--AAAAAAA
                                     "AAAAAAATCAAAAAAA",    // AAAAAAATCAAAAAAA
                                     "2=2X3=2D7=", true, true, "TC", "AATC", {12, 2, 0, 2}, {12, 2, 0, 1}},
        ExtractVariantStringTestCase{"Simple repeat masking sliding window test 6c - deletion, should not mask",
                                     "AAATCAAAAAAAAA",      // AAATCAA--AAAAAAA
                                     "AAAAAAATCAAAAAAA",    // AAAAAAATCAAAAAAA
                                     "3=2X2=2D7=", true, true, "TC", "AATC", {12, 2, 0, 2}, {12, 2, 0, 1}},
        ExtractVariantStringTestCase{"Simple repeat masking sliding window test 6d - deletion, should not mask",
                                     "AAAATCAAAAAAAA",      // AAAATCA--AAAAAAA
                                     "AAAAAAATCAAAAAAA",    // AAAAAAATCAAAAAAA
                                     "4=2X1=2D7=", true, true, "TC", "AATC", {12, 2, 0, 2}, {12, 2, 0, 1}},
        ExtractVariantStringTestCase{"Simple repeat masking sliding window test 6e - deletion, should mask",
                                     "AAAAATCAAAAAAA",      // AAAAATC--AAAAAAA
                                     "AAAAAAATCAAAAAAA",    // AAAAAAATCAAAAAAA
                                     "5=2X2D7=", true, true, "TC", "AAtc", {12, 2, 0, 0}, {12, 2, 0, 0}},
        ExtractVariantStringTestCase{"Simple repeat masking sliding window test 6f - deletion, should not mask",
                                     "AAAAAATCAAAAAA",      // AAAAAAT--CAAAAAA
                                     "AAAAAAATCAAAAAAA",    // AAAAAAATCAAAAAAA
                                     "6=1X2D1X6=", true, true, "TC", "ATCA", {12, 2, 0, 2}, {12, 2, 0, 1}},
        ExtractVariantStringTestCase{"Simple repeat masking sliding window test 6g - deletion, should mask",
                                     "AAAAAAATCAAAAA",      // AAAAAAA--TCAAAAA
                                     "AAAAAAATCAAAAAAA",    // AAAAAAATCAAAAAAA
                                     "7=2D2X5=", true, true, "TC", "tcAA", {12, 2, 0, 0}, {12, 2, 0, 0}},
        ExtractVariantStringTestCase{"Simple repeat masking sliding window test 6h - deletion, should not mask",
                                     "AAAAAAAATCAAAA",      // AAAAAAA--ATCAAAA
                                     "AAAAAAATCAAAAAAA",    // AAAAAAATCAAAAAAA
                                     "7=2D1=2X4=", true, true, "TC", "TCAA", {12, 2, 0, 2}, {12, 2, 0, 1}},
        ExtractVariantStringTestCase{"Simple repeat masking sliding window test 6i - deletion, should not mask",
                                     "AAAAAAAAATCAAA",      // AAAAAAA--AATCAAA
                                     "AAAAAAATCAAAAAAA",    // AAAAAAATCAAAAAAA
                                     "7=2D2=2X3=", true, true, "TC", "TCAA", {12, 2, 0, 2}, {12, 2, 0, 1}},
        ExtractVariantStringTestCase{"Simple repeat masking sliding window test 6j - deletion, should not mask",
                                     "AAAAAAAAAATCAA",      // AAAAAAA--AAATCAA
                                     "AAAAAAATCAAAAAAA",    // AAAAAAATCAAAAAAA
                                     "7=2D3=2X2=", true, true, "TC", "TCAA", {12, 2, 0, 2}, {12, 2, 0, 1}},
        ExtractVariantStringTestCase{"Simple repeat masking sliding window test 6k - deletion, should not mask",
                                     "AAAAAAAAAAATCA",      // AAAAAAA--AAAATCA
                                     "AAAAAAATCAAAAAAA",    // AAAAAAATCAAAAAAA
                                     "7=2D4=2X1=", true, true, "TC", "TCAA", {12, 2, 0, 2}, {12, 2, 0, 1}},
        /*
         * Test simple repeat masking in the same sequence (query).
        */
        ExtractVariantStringTestCase{"Simple repeat masking in same sequence window test 7a - should not mask",
                                     "AAAAAAATCTAAAAAA",    // AAAAAAATCTAAAAAA
                                     "AAAAAAAAAAAAAA",      // AAAAAAA--AAAAAAA
                                     "7=2I1X6=", true, true, "TCT", "A", {13, 1, 2, 0}, {13, 1, 1, 0}},
        ExtractVariantStringTestCase{"Simple repeat masking in same sequence window test 7b - should mask",
                                     "AAAAAAATCTCAAAAA",    // AAAAAAATCTCAAAAA
                                     "AAAAAAAAAAAAAA",      // AAAAAAA--AAAAAAA
                                     "7=2I2X5=", true, true, "tcTC", "AA", {12, 2, 0, 0}, {12, 2, 0, 0}},
        ExtractVariantStringTestCase{"Simple repeat masking in same sequence window test 7c - should not mask",
                                     "AAAAAATTCAAAAAAA",    // AAAAAATTCAAAAAA
                                     "AAAAAAAAAAAAAA",      // AAAAAAA--AAAAAAA
                                     "6=1X2I7=", true, true, "TTC", "A", {13, 1, 2, 0}, {13, 1, 1, 0}},
        ExtractVariantStringTestCase{"Simple repeat masking in same sequence window test 7d - should mask",
                                     "AAAAATCTCAAAAAAA",    // AAAAAAATCTCAAAAA
                                     "AAAAAAAAAAAAAA",      // AAAAAAA--AAAAAAA
                                     "5=2X2I7=", true, true, "TCtc", "AA", {12, 2, 0, 0}, {12, 2, 0, 0}},
        /*
         * Test simple repeat masking in the same sequence (target).
        */
        ExtractVariantStringTestCase{"Simple repeat masking in same sequence window test 8a - should not mask",
                                     "AAAAAAAAAAAAAA",      // AAAAAAA--AAAAAAA
                                     "AAAAAAATCTAAAAAA",    // AAAAAAATCTAAAAAA
                                     "7=2D1X6=", true, true, "A", "TCT", {13, 1, 0, 2}, {13, 1, 0, 1}},
        ExtractVariantStringTestCase{"Simple repeat masking in same sequence window test 8b - should mask",
                                     "AAAAAAAAAAAAAA",      // AAAAAAA--AAAAAAA
                                     "AAAAAAATCTCAAAAA",    // AAAAAAATCTCAAAAA
                                     "7=2D2X5=", true, true, "AA", "tcTC", {12, 2, 0, 0}, {12, 2, 0, 0}},
        ExtractVariantStringTestCase{"Simple repeat masking in same sequence window test 8c - should not mask",
                                     "AAAAAAAAAAAAAA",      // AAAAAAA--AAAAAAA
                                     "AAAAAATTCAAAAAAA",    // AAAAAATTCAAAAAA
                                     "6=1X2D7=", true, true, "A", "TTC", {13, 1, 0, 2}, {13, 1, 0, 1}},
        ExtractVariantStringTestCase{"Simple repeat masking in same sequence window test 8d - should mask",
                                     "AAAAAAAAAAAAAA",      // AAAAAAA--AAAAAAA
                                     "AAAAATCTCAAAAAAA",    // AAAAAAATCTCAAAAA
                                     "5=2X2D7=", true, true, "AA", "TCtc", {12, 2, 0, 0}, {12, 2, 0, 0}},
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
    std::vector<std::tuple<std::string, std::string, std::string, std::string, Alignment::DiffCounts, bool, bool>> testData = {
        {"Empty input", "", "", "", Alignment::DiffCounts(), false, false},
        {"Exact match", "4=", "", "", Alignment::DiffCounts(4, 0, 0, 0), false, false},
        {"Variant pos out of bounds (X)", "3=1X1=", "", "", Alignment::DiffCounts(), false, true},
        {"Variant pos out of bounds (I)", "3=1I1=", "", "", Alignment::DiffCounts(), false, true},
        {"Variant pos out of bounds (D)", "3=1D1=", "", "", Alignment::DiffCounts(), false, true},
        {"Valid case, 1X", "3=1X1=", "A", "C", Alignment::DiffCounts(4, 1, 0, 0), false, false},
        {"Valid case, 1X, masked", "3=1X1=", "a", "c", Alignment::DiffCounts(4, 0, 0, 0), false, false},
        {"Valid case, 1X, mixed mask, should throw", "3=1X1=", "a", "C", Alignment::DiffCounts(), false, true},
        {"Valid case with multiple diffs and 1 masked I", "3=1X1=1D2=2I2=1I", "ATAt", "CT", Alignment::DiffCounts(8, 1, 2, 1), false, false},
        {"Valid case with clipping and reference skip", "3S3=1X1=1D2=1I2N2=1I3H", "ATt", "CT", Alignment::DiffCounts(8, 1, 1, 1), false, false},
        {"The 2-base I has mixed masking, but the parameter says not to throw", "3=1X1=1D2=2I2=1I", "ATat", "CT", Alignment::DiffCounts(8, 1, 1, 1), false, false},
        {"The 2-base D has mixed masking, but the parameter says not to throw", "3=1X1=2D2=2I2=1I", "ATat", "CTc", Alignment::DiffCounts(8, 1, 1, 1), false, false},
        {"Bad case, the 2-base I has mixed masking", "3=1X1=1D2=2I2=1I", "ATat", "CT", Alignment::DiffCounts(), true, true},
        {"Bad case, the 2-base B has mixed masking", "3=1X1=2D2=2I2=1I", "ATat", "CTc", Alignment::DiffCounts(), true, true},
    };
    // clang-format on

    for (const auto& data : testData) {
        // Inputs.
        const std::string testName = std::get<0>(data);
        PacBio::BAM::Cigar cigar(std::get<1>(data));
        const std::string& queryVariants = std::get<2>(data);
        const std::string& targetVariants = std::get<3>(data);
        const PacBio::Pancake::Alignment::DiffCounts& expected = std::get<4>(data);
        bool throwOnPartiallyMaskedIndels = std::get<5>(data);
        bool shouldThrow = std::get<6>(data);

        // Name the test.
        SCOPED_TRACE("ComputeDiffCounts-" + testName);

        // Run.
        if (shouldThrow) {
            EXPECT_THROW(
                {
                    PacBio::Pancake::ComputeDiffCounts(cigar, queryVariants, targetVariants,
                                                       throwOnPartiallyMaskedIndels);
                },
                std::runtime_error);
        } else {
            auto result = PacBio::Pancake::ComputeDiffCounts(cigar, queryVariants, targetVariants,
                                                             throwOnPartiallyMaskedIndels);
            // Evaluate.
            EXPECT_EQ(expected, result);
        }
    }
}

TEST(Test_AlignmentTools_NormalizeAlignmentInPlace, ArrayOfTests)
{
    // clang-format off
    std::vector<std::tuple<std::string, std::string, std::string, std::string, std::string, bool>> testData = {
        {"Empty input", "", "", "", "", false},

        {"Exact match",
                        // Input alignment.
                        "ACTG",
                        "ACTG",
                        // Expected alignment.
                        "ACTG",
                        "ACTG",
                        false},

        {"Just a mismatch",
                        // Input alignment.
                        "CAC",
                        "CGC",
                        // Expected alignment.
                        "CAC",
                        "CGC",
                        false},

        {"Simple normalization with 1 indel and 1 mismatch",
                        // Input alignment.
                        "TTGACACT",
                        "TTG-TACT",
                        // Expected alignment.
                        "TTGACACT",
                        "TTGT-ACT",
                        false},
        {"Simple with two deletions in the query",
                        // Input alignment.
                        "AC--TAAC",
                        "ACTATAAC",
                        // Expected alignment.
                        "ACTA-A-C",
                        "ACTATAAC",
                        false},
        {"Test reverse complement alignment of the previous one. Shows that left alignment is not symmetric.",
                        // Input alignment.
                        "GTTATAGT",
                        "GTTA--GT",
                        // Expected alignment.
                        "GTTATAGT",
                        "GTTA--GT",
                        false},

        {"Test shifting of gaps on the query",
                        // Input alignment.
                        "-C--CGT",
                        "CCGAC-T",
                        // Expected alignment.
                        "CCG---T",
                        "CCGAC-T",                  // TODO: Clear meaningles "-":"-" alignments.
                        false},

        {"Test shifting of gaps on the target",
                        // Input alignment.
                        "ATAT-AGCCGGC",
                        "ATATTA---GGC",
                        // Expected alignment.
                        "ATAT-AGCCGGC",
                        "ATATTAG--G-C",
                        false},


    };
    // clang-format on

    for (const auto& data : testData) {
        // Inputs.
        const std::string testName = std::get<0>(data);
        const std::string& queryAln = std::get<1>(data);
        const std::string& targetAln = std::get<2>(data);
        const std::string& expectedQueryAln = std::get<3>(data);
        const std::string& expectedTargetAln = std::get<4>(data);
        const bool shouldThrow = std::get<5>(data);

        // Name the test.
        SCOPED_TRACE("NormalizeAlignmentInPlace-" + testName);

        std::string resultQueryAln = queryAln;
        std::string resultTargetAln = targetAln;

        // std::cerr << "Q: " << queryAln << "\n"
        //           << "T: " << targetAln << "\n";

        // Run.
        if (shouldThrow) {
            EXPECT_THROW(
                { PacBio::Pancake::NormalizeAlignmentInPlace(resultQueryAln, resultTargetAln); },
                std::runtime_error);
        } else {
            PacBio::Pancake::NormalizeAlignmentInPlace(resultQueryAln, resultTargetAln);
            // std::cerr << "Q: " << resultQueryAln << "\n"
            //           << "T: " << resultTargetAln << "\n";
            // Evaluate.
            EXPECT_EQ(expectedQueryAln, resultQueryAln);
            EXPECT_EQ(expectedTargetAln, resultTargetAln);
        }
        // std::cerr << "Test done.\n----------------------\n";
    }
}
}
