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
        std::vector<std::string> expectedOverlaps;
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
            {
            },
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
            {
                "000000000 000000000 -12670 99.80 0 0 12695 12695 0 602 13301 14001 contained",
                "000000000 000000001 -8210 98.58 0 4321 12695 12695 1 532 8861 9501 u",
            },
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
            {
                "000000000 000000000 -11811 100.00 0 0 11811 11811 0 0 11811 11811 contained",
                "000000000 000000001 -9053 99.74 0 0 9077 11811 0 981 10059 10059 5",
                "000000000 000000005 -7889 99.94 0 0 7895 11811 1 0 7894 11307 5",
                "000000000 000000004 -6015 99.82 0 0 6031 11811 0 2564 8590 8590 5",
                "000000000 000000002 -4942 99.68 0 0 4959 11811 0 4181 9139 9139 5",
                "000000000 000000003 -3829 99.77 0 0 3845 11811 1 0 3838 10150 5",

                "000000001 000000001 -10059 100.00 0 0 10059 10059 0 0 10059 10059 contained",
                "000000001 000000000 -9053 99.74 0 981 10059 10059 0 0 9077 11811 3",
                "000000001 000000005 -8859 99.83 0 0 8877 10059 1 0 8874 11307 5",
                "000000001 000000004 -6994 99.80 0 0 7010 10059 0 1582 8590 8590 5",
                "000000001 000000002 -5921 99.73 0 0 5937 10059 0 3201 9139 9139 5",
                "000000001 000000003 -4805 99.77 0 0 4823 10059 1 0 4816 10150 5",

                "000000002 000000002 -9139 100.00 0 0 9139 9139 0 0 9139 9139 contained",
                "000000002 000000005 -8359 99.86 0 768 9139 9139 1 2936 11307 11307 3",
                "000000002 000000003 -8010 99.90 0 0 8024 9139 1 0 8018 10150 5",
                "000000002 000000004 -7498 99.73 0 1619 9139 9139 0 0 7518 8590 3",
                "000000002 000000001 -5922 99.73 0 3201 9139 9139 0 0 5938 10059 3",
                "000000002 000000000 -4942 99.68 0 4181 9139 9139 0 0 4960 11811 3",

                "000000003 000000003 -10150 100.00 0 0 10150 10150 0 0 10150 10150 contained",
                "000000003 000000002 -8010 99.90 0 0 8018 10150 1 0 8024 9139 5",
                "000000003 000000005 -7242 99.89 0 0 7250 10150 0 4051 11307 11307 5",
                "000000003 000000004 -6384 99.78 0 0 6398 10150 1 0 6402 8590 5",
                "000000003 000000001 -4805 99.77 0 0 4816 10150 1 0 4823 10059 5",
                "000000003 000000000 -3829 99.77 0 0 3838 10150 1 0 3845 11811 5",

                "000000004 000000004 -8590 100.00 0 0 8590 8590 0 0 8590 8590 contained",
                "000000004 000000005 -8582 99.91 0 0 8590 8590 1 1864 10456 11307 contained",
                "000000004 000000002 -7497 99.73 0 0 7517 8590 0 1619 9139 9139 5",
                "000000004 000000001 -6994 99.80 0 1582 8590 8590 0 0 7012 10059 3",
                "000000004 000000003 -6384 99.78 0 0 6402 8590 1 0 6398 10150 5",
                "000000004 000000000 -6015 99.82 0 2564 8590 8590 0 0 6031 11811 3",

                "000000005 000000005 -11307 100.00 0 0 11307 11307 0 0 11307 11307 contained",
                "000000005 000000001 -8859 99.83 0 0 8874 11307 1 0 8877 10059 5",
                "000000005 000000004 -8582 99.91 0 1864 10456 11307 1 0 8590 8590 contains",
                "000000005 000000002 -8358 99.86 0 2937 11307 11307 1 768 9139 9139 3",
                "000000005 000000000 -7889 99.94 0 0 7894 11307 1 0 7895 11811 5",
                "000000005 000000003 -7242 99.89 0 4051 11307 11307 0 0 7250 10150 3",
            },
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

        SCOPED_TRACE(data.testName);

        std::cerr << "testName = " << data.testName << "\n";

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

                std::cerr << PacBio::Pancake::OverlapWriterBase::PrintOverlapAsM4(mapping, "", "",
                                                                                  true, false)
                          << "\n";

                resultsStr.emplace_back(PacBio::Pancake::OverlapWriterBase::PrintOverlapAsM4(
                    mapping, "", "", true, false));
            }
        }

        ASSERT_EQ(data.expectedOverlaps, resultsStr);
    }
}
