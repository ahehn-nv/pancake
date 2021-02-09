#include <PancakeTestData.h>
#include <gtest/gtest.h>
#include <pacbio/pancake/MapperHiFi.h>
#include <pacbio/pancake/OverlapWriterBase.h>
#include <iostream>
#include <memory>
#include "TestHelperUtils.h"

TEST(MapperHiFi, CheckMappping_LoadFromFile)
{
    struct TestData
    {
        std::string testName;
        std::string targetFile;
        std::string queryFile;
        PacBio::Pancake::SeedDB::SeedDBParameters seedParams;
        std::string expectedOverlapFile;
    };

    // clang-format off
    std::vector<TestData> testData = {
        {
            "Empty input",
            // Target.
            "",
            // Query.
            "",
            // SeedParams.
            PacBio::Pancake::SeedDB::SeedDBParameters{28, 80, 0, false, false, 255, true},
            // Expected results.
            ""
        },
        {
            "Small synthetic test case with 2 refs and one query.",
            // Target.
            PacBio::PancakeTestsConfig::Data_Dir + "/mapping/test-1-no-secondary-aln.ref.fasta",
            // Query.
            PacBio::PancakeTestsConfig::Data_Dir + "/mapping/test-1-no-secondary-aln.reads.fasta",
            // SeedParams.
            PacBio::Pancake::SeedDB::SeedDBParameters{28, 80, 0, false, false, 255, true},
            // Expected results.
            PacBio::PancakeTestsConfig::Data_Dir + "/mapping/test-1-no-secondary-aln.out.ovl"
        },
        {
            "A set of real reads, used as both the query and target sequences",
            // Target.
            PacBio::PancakeTestsConfig::Data_Dir + "/hifi-ovl/reads.pile1-5prime.fasta",
            // Query.
            PacBio::PancakeTestsConfig::Data_Dir + "/hifi-ovl/reads.pile1-5prime.fasta",
            // SeedParams.
            PacBio::Pancake::SeedDB::SeedDBParameters{28, 80, 0, false, false, 255, true},
            // Expected results.
            PacBio::PancakeTestsConfig::Data_Dir + "/hifi-ovl/reads.pile1-5prime.out.all_vs_all.ovl"
        },

    };
    // clang-format on

    PacBio::Pancake::OverlapHifiSettings settings;

    for (const auto& data : testData) {
        // Load the sequences from files, and take only the first one.
        const std::vector<std::string> allTargetSeqs =
            data.targetFile.empty()
                ? std::vector<std::string>()
                : PacBio::PancakeTests::HelperLoadFastaAsStringVector(data.targetFile);
        const std::vector<std::string> allQuerySeqs =
            data.queryFile.empty()
                ? std::vector<std::string>()
                : PacBio::PancakeTests::HelperLoadFastaAsStringVector(data.queryFile);

        std::vector<std::string> expectedOverlaps =
            data.expectedOverlapFile.empty()
                ? std::vector<std::string>()
                : PacBio::PancakeTests::HelperLoadFile(data.expectedOverlapFile);

        SCOPED_TRACE(data.testName);

        // std::cerr << "testName = " << data.testName << "\n";

        std::vector<PacBio::Pancake::OverlapHiFi::MapperResult> result =
            PacBio::Pancake::OverlapHiFi::MapHiFi(allTargetSeqs, allQuerySeqs, data.seedParams,
                                                  settings);

        std::vector<std::string> resultsStr;
        for (const auto& queryMappings : result) {
            for (size_t i = 0; i < queryMappings.overlaps.size(); ++i) {
                if (queryMappings.overlaps[i] == nullptr) {
                    continue;
                }
                const auto& mapping = queryMappings.overlaps[i];

                // std::cerr << PacBio::Pancake::OverlapWriterBase::PrintOverlapAsM4(mapping, "", "",
                //                                                                   true, false)
                //           << "\n";

                resultsStr.emplace_back(PacBio::Pancake::OverlapWriterBase::PrintOverlapAsM4(
                    mapping, "", "", true, false));
            }
        }

        ASSERT_EQ(expectedOverlaps, resultsStr);
    }
}
