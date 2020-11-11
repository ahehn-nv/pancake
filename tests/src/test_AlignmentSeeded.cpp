#include <gtest/gtest.h>

#include <pacbio/pancake/AlignmentSeeded.h>

#include <cstdint>
#include <string>
#include <tuple>
#include <vector>

using namespace PacBio::Pancake;

TEST(AlignmentSeeded, ExtractAlignmentRegions_ArrayOfTests)
{
    struct TestData
    {
        std::string testName;
        std::vector<PacBio::Pancake::SeedHit> sortedHits;
        int32_t queryLen = 0;
        int32_t targetLen = 0;
        bool isRev = false;
        int32_t minAlignmentSpan = 200;
        int32_t maxFlankExtensionDist = 5000;
        double flankExtensionFactor = 1.3;
        bool expectedThrow = false;
        std::vector<PacBio::Pancake::AlignmentRegion> expectedRegions;
    };

    // clang-format off
    std::vector<TestData> allTests{
        TestData{
            "Empty input", {}, 0, 0, false, 200, 5000, 1.3, false, {},
        },

        TestData{
            "Test parameters - invalid query length, should throw", {}, -1, 0, false, 200, 5000, 1.3, true, {},
        },
        TestData{
            "Test parameters - invalid target length, should throw", {}, 0, -1, false, 200, 5000, 1.3, true, {},
        },
        TestData{
            "Test parameters - invalid flankExtensionFactor, value is below 1.0", {}, 0, 0, false, 200, 5000, 0.8, true, {},
        },
        TestData{
            "Test parameters - invalid flankExtensionFactor, value is below 0.0", {}, 0, 0, false, 200, 5000, -1.0, true, {},
        },

        TestData{
            "Single seed - global aln is performed from seed to seed, so the final global coords start and end at the same point.",
            {   // sortedHits
                // targetId, targetRev, targetPos, queryPos, targetSpan, querySpan, flags
                SeedHit(0, false, 5, 5, 5, 5, 0),
            },
            // queryLen, targetLen, isRev, minAlignmentSpan, maxFlankExtensionDist, flankExtensionFactor
            20, 20, false, 200, 5000, 1.3,
            // expectedThrow
            false,
            // expectedRegions
            {
                {0, 5, 0, 5, false, RegionType::FRONT, 0},
                {5, 15, 5, 15, false, RegionType::BACK, 1},
            },
        },

        TestData{
            "Two seeds",
            {   // sortedHits
                // targetId, targetRev, targetPos, queryPos, targetSpan, querySpan, flags
                SeedHit(0, false, 5, 5, 5, 5, 0),
                SeedHit(0, false, 10, 10, 5, 5, 0),
            },
            // queryLen, targetLen, isRev, minAlignmentSpan, maxFlankExtensionDist, flankExtensionFactor
            20, 20, false, 200, 5000, 1.3,
            // expectedThrow
            false,
            // expectedRegions
            {
                // qStart, qSpan, tStart, tSpan, queryRev, regionType, regionId
                {0, 5, 0, 5, false, RegionType::FRONT, 0},
                {5, 5, 5, 5, false, RegionType::GLOBAL, 1},
                {10, 10, 10, 10, false, RegionType::BACK, 2},
            },
        },

        TestData{
            "Three seeds",
            {   // sortedHits
                // targetId, targetRev, targetPos, queryPos, targetSpan, querySpan, flags
                SeedHit(0, false, 5, 5, 15, 15, 0),
                SeedHit(0, false, 400, 400, 15, 15, 0),
                SeedHit(0, false, 900, 900, 15, 15, 0),
            },
            // queryLen, targetLen, isRev, minAlignmentSpan, maxFlankExtensionDist, flankExtensionFactor
            1000, 1000, false, 200, 5000, 1.3,
            // expectedThrow
            false,
            // expectedRegions
            {
                // qStart, qSpan, tStart, tSpan, queryRev, regionType, regionId
                {0, 5, 0, 5, false, RegionType::FRONT, 0},
                {5, 395, 5, 395, false, RegionType::GLOBAL, 1},
                {400, 500, 400, 500, false, RegionType::GLOBAL, 2},
                {900, 100, 900, 100, false, RegionType::BACK, 3},
            },
        },

        TestData{
            "Three seeds but middle one is skipped because distance from previous is too short",
            {   // sortedHits
                // targetId, targetRev, targetPos, queryPos, targetSpan, querySpan, flags
                SeedHit(0, false, 5, 5, 15, 15, 0),
                SeedHit(0, false, 100, 100, 15, 15, 0),
                SeedHit(0, false, 900, 900, 15, 15, 0),
            },
            // queryLen, targetLen, isRev, minAlignmentSpan, maxFlankExtensionDist, flankExtensionFactor
            1000, 1000, false, 200, 5000, 1.3,
            // expectedThrow
            false,
            // expectedRegions
            {
                // qStart, qSpan, tStart, tSpan, queryRev, regionType, regionId
                {0, 5, 0, 5, false, RegionType::FRONT, 0},
                {5, 895, 5, 895, false, RegionType::GLOBAL, 1},
                {900, 100, 900, 100, false, RegionType::BACK, 2},
            },
        },

        TestData{
            "Three seeds, same as previous, where the middle one is close to the previous hit. Here, we allow any distance for alignment.",
            {   // sortedHits
                // targetId, targetRev, targetPos, queryPos, targetSpan, querySpan, flags
                SeedHit(0, false, 5, 5, 15, 15, 0),
                SeedHit(0, false, 100, 100, 15, 15, 0),
                SeedHit(0, false, 900, 900, 15, 15, 0),
            },
            // queryLen, targetLen, isRev, minAlignmentSpan, maxFlankExtensionDist, flankExtensionFactor
            1000, 1000, false, 0, 5000, 1.3,
            // expectedThrow
            false,
            // expectedRegions
            {
                // qStart, qSpan, tStart, tSpan, queryRev, regionType, regionId
                {0, 5, 0, 5, false, RegionType::FRONT, 0},
                {5, 95, 5, 95, false, RegionType::GLOBAL, 1},
                {100, 800, 100, 800, false, RegionType::GLOBAL, 1},
                {900, 100, 900, 100, false, RegionType::BACK, 2},
            },
        },

        TestData{
            "Left flank is too far to be aligned fully, so it's limited to 5000bp",
            {   // sortedHits
                // targetId, targetRev, targetPos, queryPos, targetSpan, querySpan, flags
                SeedHit(0, false, 12000, 12000, 15, 15, 0),
                SeedHit(0, false, 14000, 14000, 15, 15, 0),
                SeedHit(0, false, 15000, 15000, 15, 15, 0),
            },
            // queryLen, targetLen, isRev, minAlignmentSpan, maxFlankExtensionDist, flankExtensionFactor
            17000, 17000, false, 200, 5000, 1.3,
            // expectedThrow
            false,
            // expectedRegions
            {
                // qStart, qSpan, tStart, tSpan, queryRev, regionType, regionId
                {7000, 5000, 7000, 5000, false, RegionType::FRONT, 0},
                {12000, 2000, 12000, 2000, false, RegionType::GLOBAL, 1},
                {14000, 1000, 14000, 1000, false, RegionType::GLOBAL, 2},
                {15000, 2000, 15000, 2000, false, RegionType::BACK, 3},
            },
        },

        TestData{
            "Right and left flanks are far from the edge, but we allow that with the maxFlankExtensionDist parameter",
            {   // sortedHits
                // targetId, targetRev, targetPos, queryPos, targetSpan, querySpan, flags
                SeedHit(0, false, 12000, 12000, 15, 15, 0),
                SeedHit(0, false, 14000, 14000, 15, 15, 0),
                SeedHit(0, false, 15000, 15000, 15, 15, 0),
            },
            // queryLen, targetLen, isRev, minAlignmentSpan, maxFlankExtensionDist, flankExtensionFactor
            25000, 25000, false, 200, -1, 1.3,
            // expectedThrow
            false,
            // expectedRegions
            {
                // qStart, qSpan, tStart, tSpan, queryRev, regionType, regionId
                {0, 12000, 0, 12000, false, RegionType::FRONT, 0},
                {12000, 2000, 12000, 2000, false, RegionType::GLOBAL, 1},
                {14000, 1000, 14000, 1000, false, RegionType::GLOBAL, 2},
                {15000, 10000, 15000, 10000, false, RegionType::BACK, 3},
            },
        },

        TestData{
            "Hits are not sorted, should throw",
            {   // sortedHits
                // targetId, targetRev, targetPos, queryPos, targetSpan, querySpan, flags
                SeedHit(0, false, 12000, 12000, 15, 15, 0),
                SeedHit(0, false, 14000, 14000, 15, 15, 0),
                SeedHit(0, false, 13000, 13000, 15, 15, 0),
                SeedHit(0, false, 15000, 15000, 15, 15, 0),
            },
            // queryLen, targetLen, isRev, minAlignmentSpan, maxFlankExtensionDist, flankExtensionFactor
            25000, 25000, false, 200, 5000, 1.3,
            // expectedThrow
            true,
            // expectedRegions
            {
            },
        },

        TestData{
            "Query is internal to target, and flank extension in target should cover 1.5x the flank length in query to allow for potential indels to be aligned",
            {   // sortedHits
                // targetId, targetRev, targetPos, queryPos, targetSpan, querySpan, flags
                SeedHit(0, false, 12000, 2000, 15, 15, 0),
                SeedHit(0, false, 14000, 4000, 15, 15, 0),
                SeedHit(0, false, 15000, 5000, 15, 15, 0),
            },
            // queryLen, targetLen, isRev, minAlignmentSpan, maxFlankExtensionDist, flankExtensionFactor
            7000, 25000, false, 200, 5000, 1.5,
            // expectedThrow
            false,
            // expectedRegions
            {
                // qStart, qSpan, tStart, tSpan, queryRev, regionType, regionId
                {0, 2000, 9000, 3000, false, RegionType::FRONT, 0},
                {2000, 2000, 12000, 2000, false, RegionType::GLOBAL, 1},
                {4000, 1000, 14000, 1000, false, RegionType::GLOBAL, 2},
                {5000, 2000, 15000, 3000, false, RegionType::BACK, 3},
            },
        },

        TestData{
            "Query is internal to target, but no extra flank extension in target is allowed",
            {   // sortedHits
                // targetId, targetRev, targetPos, queryPos, targetSpan, querySpan, flags
                SeedHit(0, false, 12000, 2000, 15, 15, 0),
                SeedHit(0, false, 14000, 4000, 15, 15, 0),
                SeedHit(0, false, 15000, 5000, 15, 15, 0),
            },
            // queryLen, targetLen, isRev, minAlignmentSpan, maxFlankExtensionDist, flankExtensionFactor
            7000, 25000, false, 200, 5000, 1.0,
            // expectedThrow
            false,
            // expectedRegions
            {
                // qStart, qSpan, tStart, tSpan, queryRev, regionType, regionId
                {0, 2000, 10000, 2000, false, RegionType::FRONT, 0},
                {2000, 2000, 12000, 2000, false, RegionType::GLOBAL, 1},
                {4000, 1000, 14000, 1000, false, RegionType::GLOBAL, 2},
                {5000, 2000, 15000, 2000, false, RegionType::BACK, 3},
            },
        },

    };
    // clang-format on

    for (const auto& data : allTests) {
        // Name the test.
        SCOPED_TRACE(data.testName);

        if (data.expectedThrow) {
            EXPECT_THROW(
                {
                    std::vector<AlignmentRegion> result = ExtractAlignmentRegions(
                        data.sortedHits, data.queryLen, data.targetLen, data.isRev,
                        data.minAlignmentSpan, data.maxFlankExtensionDist,
                        data.flankExtensionFactor);
                },
                std::runtime_error);

        } else {
            // Run the unit under test.
            std::vector<AlignmentRegion> result = ExtractAlignmentRegions(
                data.sortedHits, data.queryLen, data.targetLen, data.isRev, data.minAlignmentSpan,
                data.maxFlankExtensionDist, data.flankExtensionFactor);

            // std::cerr << "Test name: " << data.testName << "\n";
            // std::cerr << "Results:\n";
            // for (size_t i = 0; i < result.size(); ++i) {
            //     std::cerr << "[i = " << i << "] " << result[i] << "\n";
            // }
            // std::cerr << "Expected:\n";
            // for (size_t i = 0; i < data.expectedRegions.size(); ++i) {
            //     std::cerr << "[i = " << i << "] " << data.expectedRegions[i] << "\n";
            // }
            // std::cerr << "\n";

            // Evaluate.
            EXPECT_EQ(data.expectedRegions, result);
        }
    }
}

TEST(AlignmentSeeded, AlignSingleRegion_ArrayOfTests)
{
    // ASSERT_EQ(expected, result);
    // RegionsToAlignResults AlignRegionsGeneric(const RegionsToAlign& regions,
    //                                           AlignerBasePtr& alignerGlobal,
    //                                           AlignerBasePtr& alignerExt);
}

TEST(AlignmentSeeded, AlignmentSeeded_ArrayOfTests)
{
    // OverlapPtr AlignmentSeeded(const OverlapPtr& ovl, const std::vector<SeedHit>& sortedHits,
    //                            const char* targetSeq, const int32_t targetLen, const char* queryFwd,
    //                            const char* queryRev, const int32_t queryLen, int32_t minAlignmentSpan,
    //                            int32_t maxFlankExtensionDist, AlignerBasePtr& alignerGlobal,
    //                            AlignerBasePtr& alignerExt);
    // ASSERT_EQ(expected, result);
}

// Validating overlap : 000000005 000000000 105 80.50 0 -
//     137 14972 14972 1 1 15048 14902 *
//         Before : 000000005 000000000 105 0.00 0 141 14972 14972 1 0 14729 14902 *
