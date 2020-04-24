#include <gtest/gtest.h>
#include <pacbio/overlaphifi/Overlap.h>
#include <pacbio/overlaphifi/OverlapWriter.h>
#include <vector>

std::vector<std::string> WrapRunHeuristicExtendOverlapFlanks(
    const std::vector<std::string>& inOverlaps, int32_t allowedDist)
{
    std::vector<std::string> results;
    for (const auto& ovlStr : inOverlaps) {
        auto ovl = PacBio::Pancake::ParseOverlapFromString(ovlStr);
        // ovl->Type = PacBio::Pancake::DetermineOverlapType(ovl, allowedDist);
        PacBio::Pancake::HeuristicExtendOverlapFlanks(ovl, allowedDist);
        std::string resultStr =
            PacBio::Pancake::OverlapWriter::PrintOverlapAsM4(ovl, "", "", true, false);
        results.emplace_back(resultStr);
    }
    return results;
}

TEST(ExtendOverlapFlanks, EmptyInput)
{
    // Inputs.
    std::vector<PacBio::Pancake::OverlapPtr> inOverlaps = {};
    int32_t allowedDist = 50;

    // Expected results.
    std::vector<PacBio::Pancake::Overlap> expected = {};

    std::vector<PacBio::Pancake::Overlap> results;
    for (const auto& ovl : inOverlaps) {
        results.emplace_back(*PacBio::Pancake::HeuristicExtendOverlapFlanks(ovl, allowedDist));
    }

    // Evaluate.
    EXPECT_EQ(expected, results);
}

TEST(ExtendOverlapFlanks, NormalWithAllowedDist50)
{
    // Inputs.
    int32_t allowedDist = 50;
    std::vector<std::string> inOverlaps = {
        "000000001 000000002 -1000 99.00 0 0 1000 5000 0 7000 8000 8000",
        "000000001 000000002 -1000 99.00 0 0 1000 5000 1 0 1000 8000",
        "000000001 000000002 -1000 99.00 0 4000 5000 5000 1 7000 8000 8000",
        "000000001 000000002 -1000 99.00 0 4000 5000 5000 0 0 1000 8000",
        "000000001 000000002 -1000 99.00 0 4000 4998 5000 0 3 1002 5000",
        ///
        "000000001 000000002 -1000 99.00 0 4000 4998 5000 0 3 1002 5000",
        "000000001 000000002 -1000 99.00 0 5 1000 5000 0 3999 4999 5000",
        "000000001 000000002 -1000 99.00 0 3997 5000 5000 1 3998 4999 5000",
        "000000001 000000002 -1000 99.00 0 2000 3000 5000 0 10 1010 1020",
        "000000001 000000002 -1000 99.00 0 10 1010 1020 0 2000 3000 5000",
    };

    // Expected results.
    // All have 'x' types because they weren't specified in the input.
    std::vector<std::string> expected = {
        "000000001 000000002 -1000 99.00 0 0 1000 5000 0 7000 8000 8000 x",
        "000000001 000000002 -1000 99.00 0 0 1000 5000 1 0 1000 8000 x",
        "000000001 000000002 -1000 99.00 0 4000 5000 5000 1 7000 8000 8000 x",
        "000000001 000000002 -1000 99.00 0 4000 5000 5000 0 0 1000 8000 x",
        "000000001 000000002 -1000 99.00 0 3997 5000 5000 0 0 1004 5000 x",
        ///
        "000000001 000000002 -1000 99.00 0 3997 5000 5000 0 0 1004 5000 x",
        "000000001 000000002 -1000 99.00 0 0 1001 5000 0 3994 5000 5000 x",
        "000000001 000000002 -1000 99.00 0 3996 5000 5000 1 3998 5000 5000 x",
        "000000001 000000002 -1000 99.00 0 1990 3010 5000 0 0 1020 1020 x",
        "000000001 000000002 -1000 99.00 0 0 1020 1020 0 1990 3010 5000 x",
    };

    // Run and collect results.
    std::vector<std::string> results = WrapRunHeuristicExtendOverlapFlanks(inOverlaps, allowedDist);

    // Evaluate.
    EXPECT_EQ(expected, results);
}

TEST(ExtendOverlapFlanks, NormalWithAllowedDist0)
{
    /*
     * The output should be the same as input because the allowedDist is 0.
    */
    // Inputs.
    int32_t allowedDist = 0;
    std::vector<std::string> inOverlaps = {
        "000000001 000000002 -1000 99.00 0 4000 4998 5000 0 3 1002 5000",
        "000000001 000000002 -1000 99.00 0 5 1000 5000 0 3999 4999 5000",
        "000000001 000000002 -1000 99.00 0 3997 5000 5000 1 3998 4999 5000",
        "000000001 000000002 -1000 99.00 0 2000 3000 5000 0 10 1010 1020",
    };

    // Expected results.
    // All have 'x' types because they weren't specified in the input.
    std::vector<std::string> expected = {
        "000000001 000000002 -1000 99.00 0 4000 4998 5000 0 3 1002 5000 x",
        "000000001 000000002 -1000 99.00 0 5 1000 5000 0 3999 4999 5000 x",
        "000000001 000000002 -1000 99.00 0 3997 5000 5000 1 3998 4999 5000 x",
        "000000001 000000002 -1000 99.00 0 2000 3000 5000 0 10 1010 1020 x",
    };

    // Run and collect results.
    std::vector<std::string> results = WrapRunHeuristicExtendOverlapFlanks(inOverlaps, allowedDist);

    // Evaluate.
    EXPECT_EQ(expected, results);
}
