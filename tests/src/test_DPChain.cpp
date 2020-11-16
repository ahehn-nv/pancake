#include <gtest/gtest.h>

#include <pacbio/pancake/DPChain.h>
#include <cstdint>
#include <iostream>
#include <vector>

using namespace PacBio::Pancake;

TEST(DPChain, GroupByTargetAndStrand_ArrayOfTests)
{
    struct TestData
    {
        std::string testName;
        std::vector<PacBio::Pancake::SeedHit> sortedHits;
        bool expectedThrow = false;
        std::vector<PacBio::Pancake::Range> expectedResult;
    };

    // clang-format off
    std::vector<TestData> allTests{
        TestData{
            "Empty input", {}, false, {}
        },

        TestData{
            "Single point",
            {
                SeedHit(0, false, 1000, 1000, 15, 15, 0),
            },
            false,
            {
                {0, 1},
            }
        },

        TestData{
            "Two points, different target",
            {
                SeedHit(0, false, 1000, 1000, 15, 15, 0),
                SeedHit(1, false, 2000, 2000, 15, 15, 0),
            },
            false,
            {
                {0, 1},
                {1, 2},
            }
        },

        TestData{
            "Two points, same target, different strand",
            {
                SeedHit(0, false, 1000, 1000, 15, 15, 0),
                SeedHit(0, true, 2000, 2000, 15, 15, 0),
            },
            false,
            {
                {0, 1},
                {1, 2},
            }
        },

        TestData{
            "Multiple points, same target and strand",
            {
                SeedHit(0, false, 1000, 1000, 15, 15, 0),
                SeedHit(0, false, 2000, 2000, 15, 15, 0),
                SeedHit(0, false, 3000, 3000, 15, 15, 0),
                SeedHit(0, false, 4000, 4000, 15, 15, 0),
                SeedHit(0, false, 5000, 5000, 15, 15, 0),
            },
            false,
            {
                {0, 5},
            }
        },

        TestData{
            "Multiple points, three targets and strands",
            {
                SeedHit(0, false, 1000, 1000, 15, 15, 0),
                SeedHit(0, false, 2000, 2000, 15, 15, 0),
                SeedHit(0, false, 3000, 3000, 15, 15, 0),
                SeedHit(1, false, 4000, 4000, 15, 15, 0),
                SeedHit(1, false, 5000, 5000, 15, 15, 0),
                SeedHit(1, true, 6000, 6000, 15, 15, 0),
                SeedHit(2, true, 7000, 7000, 15, 15, 0),
                SeedHit(2, true, 8000, 8000, 15, 15, 0),
                SeedHit(3, true, 9000, 9000, 15, 15, 0),
            },
            false,
            {
                {0, 3},
                {3, 5},
                {5, 6},
                {6, 8},
                {8, 9},
            }
        },

    };
    // clang-format on

    for (const auto& data : allTests) {
        // Name the test.
        SCOPED_TRACE(data.testName);

        if (data.expectedThrow) {
            EXPECT_THROW({ GroupByTargetAndStrand(data.sortedHits); }, std::runtime_error);

        } else {
            // Run the unit under test.
            std::vector<PacBio::Pancake::Range> result = GroupByTargetAndStrand(data.sortedHits);

            // Evaluate.
            EXPECT_EQ(data.expectedResult, result);
        }
    }
}

TEST(DPChain, DiagonalGroup_ArrayOfTests)
{
    struct TestData
    {
        std::string testName;
        std::vector<PacBio::Pancake::SeedHit> sortedHits;
        int32_t chainBandwidth = 0;
        bool overlappingWindows = 0;
        bool expectedThrow = false;
        std::vector<PacBio::Pancake::Range> expectedResult;
    };

    // clang-format off
    std::vector<TestData> allTests{
        TestData{
            "Empty input", {}, 500, true, false, {}
        },

        TestData{
            "Single point",
            // sortedHits
            {
                SeedHit(0, false, 1000, 1000, 15, 15, 0),
            },
            // chainBandwidth, overlappingWindows
            500, false,
            // expectedThrow
            false,
            // expectedResult
            {
                {0, 1},
            }
        },

        TestData{
            "Multiple points, same target and strand",
            // sortedHits
            {
                SeedHit(0, false, 1000, 1000, 15, 15, 0),
                SeedHit(0, false, 2000, 2000, 15, 15, 0),
                SeedHit(0, false, 3000, 3000, 15, 15, 0),
                SeedHit(0, false, 4000, 4000, 15, 15, 0),
                SeedHit(0, false, 5000, 5000, 15, 15, 0),
            },
            // chainBandwidth, overlappingWindows
            500, false,
            // expectedThrow
            false,
            // expectedResult
            {
                {0, 5},
            }
        },

        TestData{
            "Multiple points, three targets and strands",
            // sortedHits
            {
                SeedHit(0, false, 1000, 1000, 15, 15, 0),
                SeedHit(0, false, 2000, 2000, 15, 15, 0),
                SeedHit(0, false, 3000, 3000, 15, 15, 0),
                SeedHit(1, false, 4000, 4000, 15, 15, 0),
                SeedHit(1, false, 5000, 5000, 15, 15, 0),
                SeedHit(1, true, 6000, 6000, 15, 15, 0),
                SeedHit(2, true, 7000, 7000, 15, 15, 0),
                SeedHit(2, true, 8000, 8000, 15, 15, 0),
                SeedHit(3, true, 9000, 9000, 15, 15, 0),
            },
            // chainBandwidth, overlappingWindows
            500, false,
            // expectedThrow
            false,
            // expectedResult
            {
                {0, 3},
                {3, 5},
                {5, 6},
                {6, 8},
                {8, 9},
            }
        },

        TestData{
            "Multiple points, same target and strand, chain bandwidth broken",
            // sortedHits
            {
                SeedHit(0, false, 1000, 1000, 15, 15, 0),
                SeedHit(0, false, 2000, 2000, 15, 15, 0),
                SeedHit(0, false, 3000, 3000, 15, 15, 0),
                SeedHit(0, false, 4700, 4000, 15, 15, 0),
                SeedHit(0, false, 5700, 5000, 15, 15, 0),
            },
            // chainBandwidth, overlappingWindows
            500, false,
            // expectedThrow
            false,
            // expectedResult
            {
                {0, 3},
                {3, 5},
            }
        },

        TestData{
            "Multiple points, multiple targets and strands, chain bandwidths broken",
            // sortedHits
            {
                SeedHit(0, false, 1000, 1000, 15, 15, 0),
                SeedHit(0, false, 2000, 2000, 15, 15, 0),
                SeedHit(0, false, 3000, 3000, 15, 15, 0),
                SeedHit(0, false, 4700, 4000, 15, 15, 0),
                SeedHit(0, false, 5700, 5000, 15, 15, 0),

                SeedHit(1, false, 1000, 1000, 15, 15, 0),
                SeedHit(1, false, 2000, 2000, 15, 15, 0),
                SeedHit(1, false, 3000, 3000, 15, 15, 0),
                SeedHit(1, false, 4700, 4000, 15, 15, 0),
                SeedHit(1, false, 5700, 5700, 15, 15, 0),
            },
            // chainBandwidth, overlappingWindows
            500, false,
            // expectedThrow
            false,
            // expectedResult
            {
                {0, 3},
                {3, 5},
                {5, 8},
                {8, 9},
                {9, 10},
            }
        },

        TestData{
            "Overlapping windows. Some hits are close to both the first and second range, and are included in both.",
            // sortedHits
            {
                // First diagonal range: from 0bp to 400bp.
                SeedHit(0, false, 1000, 1000, 15, 15, 0),
                SeedHit(0, false, 2000, 2000, 15, 15, 0),
                SeedHit(0, false, 3000, 3000, 15, 15, 0),
                SeedHit(0, false, 4000, 4000, 15, 15, 0),
                // These two seed hits should be included in the first
                // and the second range.
                SeedHit(0, false, 4300, 4000, 15, 15, 0),
                SeedHit(0, false, 4400, 4000, 15, 15, 0),

                // Second diagonal range: 700bp.
                SeedHit(0, false, 4700, 4000, 15, 15, 0),
                SeedHit(0, false, 5700, 5000, 15, 15, 0),
                SeedHit(0, false, 6700, 6000, 15, 15, 0),

                // Third diagonal range: 0bp again.
                SeedHit(0, false, 7700, 7700, 15, 15, 0),
                SeedHit(0, false, 8700, 8700, 15, 15, 0),
                SeedHit(0, false, 9700, 9700, 15, 15, 0),
            },
            // chainBandwidth, overlappingWindows
            500, true,
            // expectedThrow
            false,
            // expectedResult
            {
                {0, 6},
                {4, 9},
                {9, 12},
            }
        },
    };
    // clang-format on

    for (const auto& data : allTests) {
        // Name the test.
        SCOPED_TRACE(data.testName);

        if (data.expectedThrow) {
            EXPECT_THROW(
                { DiagonalGroup(data.sortedHits, data.chainBandwidth, data.overlappingWindows); },
                std::runtime_error);

        } else {
            // Run the unit under test.
            std::vector<PacBio::Pancake::Range> result =
                DiagonalGroup(data.sortedHits, data.chainBandwidth, data.overlappingWindows);

            // Evaluate.
            EXPECT_EQ(data.expectedResult, result);
        }
    }
}
