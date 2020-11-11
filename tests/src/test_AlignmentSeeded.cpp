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
        int32_t expectedGlobalQueryStart = 0;
        int32_t expectedGlobalQueryEnd = 0;
        int32_t expectedGlobalTargetStart = 0;
        int32_t expectedGlobalTargetEnd = 0;
    };

    // clang-format off
    std::vector<TestData> allTests{
        TestData{
            "Empty input", {}, 0, 0, false, 200, 5000, 1.3, false, {}, 0, 0, 0, 0,
        },

        TestData{
            "Single seed - global aln is performed from seed to seed, so the final global coords start and end at the same point.",
            {   // sortedHits
                // targetId, targetRev, targetPos, queryPos, targetSpan, querySpan, flags
                SeedHit(0, false, 5, 5, 5, 5, 0),
            },
            20,                             // queryLen
            20,                             // targetLen
            false,                          // isRev
            200,                            // minAlignmentSpan
            5000,                           // maxFlankExtensionDist
            1.3,                            // flankExtensionFactor
            false,                          // expectedThrow
            {   // AlignmentRegion objects
                {0, 5, 0, 5, false, RegionType::FRONT, 0},
                {5, 15, 5, 15, false, RegionType::BACK, 1},
            },
            // globalAlnQueryStart, globalAlnQueryEnd, globalAlnTargetStart, globalAlnTargetEnd
            5, 5, 5, 5,
        },

        TestData{
            "Two seeds",
            {   // sortedHits
                // targetId, targetRev, targetPos, queryPos, targetSpan, querySpan, flags
                SeedHit(0, false, 5, 5, 5, 5, 0),
                SeedHit(0, false, 10, 10, 5, 5, 0),
            },
            20,                             // queryLen
            20,                             // targetLen
            false,                          // isRev
            200,                            // minAlignmentSpan
            5000,                           // maxFlankExtensionDist
            1.3,                            // flankExtensionFactor
            false,                          // expectedThrow
            {   // AlignmentRegion objects
                // qStart, qSpan, tStart, tSpan, queryRev, regionType, regionId
                {0, 5, 0, 5, false, RegionType::FRONT, 0},
                {5, 5, 5, 5, false, RegionType::GLOBAL, 1},
                {10, 10, 10, 10, false, RegionType::BACK, 2},
            },
            // globalAlnQueryStart, globalAlnQueryEnd, globalAlnTargetStart, globalAlnTargetEnd
            5, 10, 5, 10,
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

            std::cerr << "Test name: " << data.testName << "\n";
            std::cerr << "Results:\n";
            for (size_t i = 0; i < result.size(); ++i) {
                std::cerr << "[i = " << i << "] " << result[i] << "\n";
            }

            std::cerr << "Expected:\n";
            for (size_t i = 0; i < data.expectedRegions.size(); ++i) {
                std::cerr << "[i = " << i << "] " << data.expectedRegions[i] << "\n";
            }
            // std::cerr << "data.expectedGlobalQueryStart = " << data.expectedGlobalQueryStart
            //           << "\n";
            // std::cerr << "data.expectedGlobalQueryEnd = " << data.expectedGlobalQueryEnd << "\n";
            // std::cerr << "data.expectedGlobalTargetStart = " << data.expectedGlobalTargetStart
            //           << "\n";
            // std::cerr << "data.expectedGlobalTargetEnd = " << data.expectedGlobalTargetEnd << "\n";

            // std::cerr << "result.globalAlnQueryStart = " << result.globalAlnQueryStart << "\n";
            // std::cerr << "result.globalAlnQueryEnd = " << result.globalAlnQueryEnd << "\n";
            // std::cerr << "result.globalAlnTargetStart = " << result.globalAlnTargetStart << "\n";
            // std::cerr << "result.globalAlnTargetEnd = " << result.globalAlnTargetEnd << "\n";
            std::cerr << "\n";

            // Evaluate.
            EXPECT_EQ(data.expectedRegions, result);
            // EXPECT_EQ(data.expectedGlobalQueryStart, result.globalAlnQueryStart);
            // EXPECT_EQ(data.expectedGlobalQueryEnd, result.globalAlnQueryEnd);
            // EXPECT_EQ(data.expectedGlobalTargetStart, result.globalAlnTargetStart);
            // EXPECT_EQ(data.expectedGlobalTargetEnd, result.globalAlnTargetEnd);
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
