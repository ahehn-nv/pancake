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
}
