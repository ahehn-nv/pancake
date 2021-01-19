#include <PancakeTestData.h>
#include <gtest/gtest.h>
#include <pacbio/pancake/AlignerBatchGPU.h>
#include <iostream>
#include "TestHelperUtils.h"

TEST(AlignerBatchGPU, ArrayOfTests_Small)
{
    using namespace PacBio::Pancake;

    struct TestData
    {
        std::string testName;
        std::vector<std::pair<std::string, std::string>> batchData;
        double maxMemoryFraction;
        int64_t maxMemoryCap;
        std::vector<PacBio::Pancake::AlignmentResult> expectedAlns;
    };

    // clang-format off
    std::vector<TestData> testData = {
        {
            "Batch of multiple query/target pairs.",
            {
                {"ACTGACTGAC", "ACTGTCTGAC"},
                {"ACTG", "ACTG"},
                {"A", "T"},
            },
            // Maximum memory fraction of total free memory.
            0.50,
            // Maximum memory cap.
            1024 * 1024,
            // Expected results.
            {
                PacBio::Pancake::AlignmentResult{PacBio::BAM::Cigar("4=1X5="), 10, 10, 10, 10, true, 14, 14, false},
                PacBio::Pancake::AlignmentResult{PacBio::BAM::Cigar("4="), 4, 4, 4, 4, true, 8, 8, false},
                PacBio::Pancake::AlignmentResult{PacBio::BAM::Cigar("1X"), 1, 1, 1, 1, true, -4, -4, false},
            },
        },
    };
    // clang-format on

    uint32_t maxBandwidth = 2000;
    uint32_t deviceId = 0;
    PacBio::Pancake::AlignmentParameters alnParams;

    for (const auto& data : testData) {
        // Debug info.
        SCOPED_TRACE(data.testName);
        std::cerr << "testName = " << data.testName << "\n";

        auto aligner = PacBio::Pancake::AlignerBatchGPU(alnParams, maxBandwidth, deviceId,
                                                        data.maxMemoryFraction, data.maxMemoryCap);

        for (const auto& seqPair : data.batchData) {
            const auto& query = seqPair.first;
            const auto& target = seqPair.second;
            aligner.AddSequencePair(query.c_str(), query.size(), target.c_str(), target.size());
        }

        // Run alignment.
        aligner.AlignAll();

        const std::vector<PacBio::Pancake::AlignmentResult>& results = aligner.GetAlnResults();

        std::cerr << "results.size() = " << results.size() << "\n";
        for (size_t i = 0; i < results.size(); ++i) {
            const auto& aln = results[i];
            std::cerr << "[result " << i << "] " << aln << "\n";
        }

        // Evaluate.
        ASSERT_EQ(data.expectedAlns, results);
    }
}
