#include <PancakeTestData.h>
#include <gtest/gtest.h>
#include <pacbio/pancake/MapperCLR.h>
#include <pacbio/pancake/OverlapWriterBase.h>
#include <iostream>
#include <memory>
#include "TestHelperUtils.h"

TEST(MapperCLR, CheckMappping_LoadFromFile)
{
    struct TestData
    {
        std::string testName;
        std::string targetFile;
        std::string queryFile;
        PacBio::Pancake::SeedDB::SeedDBParameters seedParamsPrimary;
        PacBio::Pancake::SeedDB::SeedDBParameters seedParamsFallback;
        std::vector<std::string> expectedOverlaps;
    };

    // clang-format off
    std::vector<TestData> testData = {
        {
            "Pair of real subreads. CCS.",
            // Target.
            PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-2-real-insert-target.fasta",
            // Query.
            PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-2-real-insert-query.fasta",
            // SeedParams - primary.
            PacBio::Pancake::SeedDB::SeedDBParameters{19, 10, 0, false, false, 255, true},
            // SeedParams - fallback.
            PacBio::Pancake::SeedDB::SeedDBParameters{10, 5, 0, false, false, 255, true},
            // Expected results.
            {
                "000000000 000000000 16398 81.12 0 0 18779 18779 0 0 18864 18865 *"
            },
        },
        {
            "RealHifiReads, subsequences with a full global match.",
            // Target.
            PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-3-real-subseq-full-global-match-target.fasta",
            // Query.
            PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-3-real-subseq-full-global-match-query.fasta",
            // SeedParams - primary.
            PacBio::Pancake::SeedDB::SeedDBParameters{19, 10, 0, false, false, 255, true},
            // SeedParams - fallback.
            PacBio::Pancake::SeedDB::SeedDBParameters{10, 5, 0, false, false, 255, true},
            // Expected results.
            {
                "000000000 000000000 15684 98.71 0 0 8111 8111 0 0 8138 8138 *"
            },
        },
        {
            "Small reverse complement perfect alignment",
            // Target.
            PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-4-small-revcmp-perfect-aln-target.fasta",
            // Query.
            PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-4-small-revcmp-perfect-aln-query.fasta",
            // SeedParams - primary.
            PacBio::Pancake::SeedDB::SeedDBParameters{19, 10, 0, false, true, 255, true},
            // SeedParams - fallback.
            PacBio::Pancake::SeedDB::SeedDBParameters{10, 5, 0, false, true, 255, true},
            // Expected results.
            {
                "000000000 000000000 300 100.00 0 0 150 150 1 0 150 150 *"
            },
        },
        {
            "Reverse complement alignment. This one was challenging for the AlignmentSeeded implementation, because offseting the reverse coordinates had a bug.",
            // Target.
            PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-5-revcmp-aln-target.fasta",
            // Query.
            PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-5-revcmp-aln-query.fasta",
            // SeedParams - primary.
            PacBio::Pancake::SeedDB::SeedDBParameters{19, 10, 0, false, false, 255, true},
            // SeedParams - fallback.
            PacBio::Pancake::SeedDB::SeedDBParameters{10, 5, 0, false, false, 255, true},
            // Expected results.
            {
                "000000000 000000000 12740 74.58 0 0 21382 22015 1 701 22028 22028 *"
            },
        },
        {
            "Small reverse complement imperfect alignment. There are 5 insertions sprinkled around, all 'A' bases.",
            // Target.
            PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-6-small-revcmp-imperfect-aln-target.fasta",
            // Query.
            PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-6-small-revcmp-imperfect-aln-query.fasta",
            // SeedParams - primary.
            PacBio::Pancake::SeedDB::SeedDBParameters{19, 10, 0, false, true, 255, true},
            // SeedParams - fallback.
            PacBio::Pancake::SeedDB::SeedDBParameters{10, 5, 0, false, true, 255, true},
            // Expected results.
            {
                "000000000 000000000 270 96.77 0 0 155 155 1 0 150 150 *"
            },
        },
        {
            "Test seed fallback. Without fallback, using large parameters (k = 19, w = 10, spacing = 0, HPC = true) will not map.",
            // Target.
            PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-7-seed-fallback-target.fasta",
            // Query.
            PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-7-seed-fallback-query.fasta",
            // SeedParams - primary.
            PacBio::Pancake::SeedDB::SeedDBParameters{19, 10, 0, false, true, 255, true},
            // SeedParams - fallback.
            PacBio::Pancake::SeedDB::SeedDBParameters{10, 5, 0, false, true, 255, true},
            // Expected results.
            {
                "000000000 000000000 392 71.07 0 0 1241 1241 0 0 1225 1225 *"
            },
        },
        {
            "The same test but without seed fallback. Without fallback, using large parameters (k = 19, w = 10, spacing = 0, HPC = true) will not map.",
            // Target.
            PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-7-seed-fallback-target.fasta",
            // Query.
            PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-7-seed-fallback-query.fasta",
            // SeedParams - primary.
            PacBio::Pancake::SeedDB::SeedDBParameters{19, 10, 0, false, true, 255, true},
            // SeedParams - fallback - same as primary.
            PacBio::Pancake::SeedDB::SeedDBParameters{19, 10, 0, false, true, 255, true},
            // Expected results.
            {
            },
        },
        {
            "RealHifiReads, poor matching pair of reads. CCS.",
            // Target.
            PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-1-poor-aln-overlapping-seeds.target.fasta",
            // Query.
            PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-1-poor-aln-overlapping-seeds.query.fasta",
            // SeedParams - primary.
            PacBio::Pancake::SeedDB::SeedDBParameters{19, 10, 0, false, false, 255, true},
            // SeedParams - fallback.
            PacBio::Pancake::SeedDB::SeedDBParameters{10, 5, 0, false, false, 255, true},
            // Expected results.
            {
                "000000000 000000000 284 68.58 0 6209 7938 43446 0 7261 8999 46238 *"
            },
        },

    };
    // clang-format on

    PacBio::Pancake::MapperCLRSettings settings;
    settings.map.freqPercentile = 0.000;

    for (const auto& data : testData) {
        // Load the sequences from files, and take only the first one.
        const std::vector<PacBio::BAM::FastaSequence> allTargetSeqs =
            PacBio::PancakeTests::HelperLoadFasta(data.targetFile);
        const std::vector<PacBio::BAM::FastaSequence> allQuerySeqs =
            PacBio::PancakeTests::HelperLoadFasta(data.queryFile);
        const std::string& target = allTargetSeqs[0].Bases();
        const std::string& query = allQuerySeqs[0].Bases();

        SCOPED_TRACE(data.testName);

        std::cerr << "testName = " << data.testName << "\n";

        settings.map.seedParams = data.seedParamsPrimary;
        settings.map.seedParamsFallback = data.seedParamsFallback;
        PacBio::Pancake::MapperCLR mapper(settings);

        std::vector<PacBio::Pancake::MapperBaseResult> result =
            mapper.MapAndAlign({target}, {query});
        std::vector<std::string> resultsStr;
        for (const auto& queryMappings : result) {
            for (const auto& mapping : queryMappings.mappings) {
                std::cerr << PacBio::Pancake::OverlapWriterBase::PrintOverlapAsM4(
                                 *mapping->mapping, "", "", true, false)
                          << "\n";

                resultsStr.emplace_back(PacBio::Pancake::OverlapWriterBase::PrintOverlapAsM4(
                    *mapping->mapping, "", "", true, false));
            }
        }

        EXPECT_EQ(data.expectedOverlaps, resultsStr);
    }
}

void DebugPrintChainedRegion(std::ostream& oss, int32_t regionId,
                             const PacBio::Pancake::ChainedRegion& cr)
{
    oss << "[regionId " << regionId << "] chain.hits = " << cr.chain.hits.size()
        << ", chain.score = " << cr.chain.score << ", chain.covQ = " << cr.chain.coveredBasesQuery
        << ", chain.covT = " << cr.chain.coveredBasesTarget << ", priority = " << cr.priority
        << ", isSuppl = " << (cr.isSupplementary ? "true" : "false") << ", ovl: "
        << PacBio::Pancake::OverlapWriterBase::PrintOverlapAsM4(*cr.mapping, "", "", true, false)
        << ", diagStart = " << (cr.mapping->Astart - cr.mapping->Bstart)
        << ", diagEnd = " << (cr.mapping->Aend - cr.mapping->Bend);

    oss << "\n";
    for (size_t j = 0; j < cr.chain.hits.size(); ++j) {
        const auto& hit = cr.chain.hits[j];
        std::cerr << "    [hit " << j << "] " << hit << "\n";
    }
}

TEST(MapperCLR, CheckMappingAndSeedHits)
{
    using namespace PacBio::Pancake;

    struct TestData
    {
        std::string testName;
        std::string targetFile;
        std::string queryFile;
        PacBio::Pancake::SeedDB::SeedDBParameters seedParamsPrimary;
        PacBio::Pancake::SeedDB::SeedDBParameters seedParamsFallback;
        std::vector<std::vector<SeedHit>> expectedSeedHits;
        std::vector<std::string> expectedMappings;
    };

    // clang-format off
    std::vector<TestData> testData = {
        {
            "RealHifiReads_Subsequences",
            // Target
            PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-1-poor-aln-overlapping-seeds.target.fasta",
            // Query.
            PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-1-poor-aln-overlapping-seeds.query.fasta",
            // SeedParams - primary
            PacBio::Pancake::SeedDB::SeedDBParameters{19, 10, 0, false, true, 100, true},
            // SeedParams - fallback
            PacBio::Pancake::SeedDB::SeedDBParameters{15, 5, 0, false, true, 100, true},
            ///// These were the old 'wrong' mappings, where we had two seed hits which had the
            ///// same X-coordinate.
            // // Expected seed hits.
            // {
            //     {
            //         SeedHit(0, false, 38493, 6585, 30, 25, 0),
            //         SeedHit(0, false, 40327, 6585, 24, 25, 0),
            //     },
            //     {
            //         SeedHit(0, false, 11779, 29565, 28, 29, 0),
            //         SeedHit(0, false, 11791, 29577, 25, 25, 0),
            //         SeedHit(0, false, 11849, 29637, 24, 27, 0),
            //     },
            // },
            // // Expected mappings.
            // {
            //     "000000000 000000000 -862 25.72 0 6468 7589 43446 0 38370 41325 46238 *",
            //     "000000000 000000000 -183 73.89 0 29527 29753 43446 0 11747 11966 46238 *"
            // },
            /////
            ///// This is what the seed hits and mappings look like after the LongMergeChains_ strictly
            ///// disallowed the overlap, even at the same position.
            // Expected seed hits.
            {
                {
                    SeedHit(0, false, 7919, 6832, 28, 28, 0),
                    SeedHit(0, false, 11021, 9922, 29, 30, 0),
                    SeedHit(0, false, 11401, 10312, 25, 24, 0),
                },
                {
                    SeedHit(0, false, 37479, 5606, 26, 24, 0),
                    SeedHit(0, false, 37486, 5612, 26, 25, 0),
                },
            },
            // Expected mappings.
            {
                "000000000 000000000 568 67.66 0 6209 12358 43446 0 7261 13474 46238 *",
                "000000000 000000000 188 71.66 0 5287 5940 43446 0 37146 37820 46238 *",
            },


        },
    };
    // clang-format on

    PacBio::Pancake::MapperCLRSettings settings;
    settings.map.freqPercentile = 0.000;

    for (const auto& data : testData) {
        // Load the sequences from files, and take only the first one.
        const std::vector<PacBio::BAM::FastaSequence> allTargetSeqs =
            PacBio::PancakeTests::HelperLoadFasta(data.targetFile);
        const std::vector<PacBio::BAM::FastaSequence> allQuerySeqs =
            PacBio::PancakeTests::HelperLoadFasta(data.queryFile);
        const std::string& target = allTargetSeqs[0].Bases();
        const std::string& query = allQuerySeqs[0].Bases();

        SCOPED_TRACE(data.testName);

        std::cerr << "testName = " << data.testName << "\n";

        settings.map.seedParams = data.seedParamsPrimary;
        settings.map.seedParamsFallback = data.seedParamsFallback;
        // settings.align = false;
        PacBio::Pancake::MapperCLR mapper(settings);

        std::vector<PacBio::Pancake::MapperBaseResult> result =
            mapper.MapAndAlign({target}, {query});
        std::vector<std::vector<PacBio::Pancake::SeedHit>> resultSeedHits;
        std::vector<std::string> resultsMappings;
        for (const auto& queryMappings : result) {
            for (const auto& mapping : queryMappings.mappings) {
                std::cerr << PacBio::Pancake::OverlapWriterBase::PrintOverlapAsM4(
                                 *mapping->mapping, "", "", true, false)
                          << "\n";

                // DebugPrintChainedRegion(std::cerr, 0, *mapping);
                resultSeedHits.emplace_back(mapping->chain.hits);

                resultsMappings.emplace_back(PacBio::Pancake::OverlapWriterBase::PrintOverlapAsM4(
                    *mapping->mapping, "", "", true, false));
            }
        }

        EXPECT_EQ(data.expectedSeedHits, resultSeedHits);
        EXPECT_EQ(data.expectedMappings, resultsMappings);
    }

    // exit(1);
}

TEST(MapperCLR, CondenseMappings)
{
    using namespace PacBio::Pancake;

    struct TestData
    {
        std::string testName;
        // tuple: <Overlap, priority, isSupplementary, isOverlapNullptr>
        std::vector<std::tuple<std::string, int32_t, bool, bool>> inMappings;
        int32_t bestNSecondary = 0;
        std::vector<std::tuple<std::string, int32_t, bool, bool>> expected;
    };

    // clang-format off
    std::vector<TestData> testData = {
        {
            "Empty input",
            {
            },
            0,
            {
            },
        },
        {
            "Single overlap, keeper",
            {
                {"000000000 000000001 -1000 99.90 0 5000 10000 10000 0 0 5000 10000 u", 0, false, false},
            },
            0,
            {
                {"000000000 000000001 -1000 99.90 0 5000 10000 10000 0 0 5000 10000 u", 0, false, false},
            },
        },
        {
            "Single overlap, nullptr",
            {
                {"", 0, false, true},
            },
            0,
            {
            },
        },
        {
            "Multiple nullptr overlap",
            {
                {"", 0, false, true},
                {"", 0, false, true},
                {"", 0, false, true},
            },
            0,
            {
            },
        },
        {
            "Multiple keeper overlaps",
            {
                {"000000000 000000001 -1000 99.90 0 5000 10000 10000 0 0 5000 10000 u", 0, false, false},
                {"000000000 000000002 -1000 99.90 0 5000 10000 10000 0 0 5000 10000 u", 0, false, false},
                {"000000000 000000003 -1000 99.90 0 5000 10000 10000 0 0 5000 10000 u", 0, false, false},
            },
            0,
            {
                {"000000000 000000001 -1000 99.90 0 5000 10000 10000 0 0 5000 10000 u", 0, false, false},
                {"000000000 000000002 -1000 99.90 0 5000 10000 10000 0 0 5000 10000 u", 0, false, false},
                {"000000000 000000003 -1000 99.90 0 5000 10000 10000 0 0 5000 10000 u", 0, false, false},
            },
        },
        {
            "One nullptr overlap to remove",
            {
                {"000000000 000000001 -1000 99.90 0 5000 10000 10000 0 0 5000 10000 u", 0, false, false},
                {"", 0, false, true},
                {"000000000 000000002 -1000 99.90 0 5000 10000 10000 0 0 5000 10000 u", 0, false, false},
            },
            0,
            {
                {"000000000 000000001 -1000 99.90 0 5000 10000 10000 0 0 5000 10000 u", 0, false, false},
                {"000000000 000000002 -1000 99.90 0 5000 10000 10000 0 0 5000 10000 u", 0, false, false},
            },
        },
        {
            "Filtering all secondary overlaps",
            {
                {"000000000 000000001 -1000 99.90 0 5000 10000 12000 0 0 5000 10000 u", 0, false, false},       // Primary
                {"000000000 000000001 -1000 99.90 0 10000 12000 12000 0 7000 9000 10000 u", 0, true, false},    // Supplementary
                {"000000000 000000002 -1000 99.90 0 5000 10000 10000 0 0 5000 10000 u", 1, false, false},       // Secondary
                {"000000000 000000003 -1000 99.90 0 5000 10000 10000 0 0 5000 10000 u", 1, false, false},       // Secondary
                {"000000000 000000004 -1000 99.90 0 5000 10000 10000 0 0 5000 10000 u", 1, false, false},       // Secondary
                {"000000005 000000006 -1000 99.90 0 5000 10000 12000 0 0 5000 10000 u", 0, false, false},       // Another primary
            },
            0,
            {
                {"000000000 000000001 -1000 99.90 0 5000 10000 12000 0 0 5000 10000 u", 0, false, false},       // Primary
                {"000000000 000000001 -1000 99.90 0 10000 12000 12000 0 7000 9000 10000 u", 0, true, false},    // Supplementary
                {"000000005 000000006 -1000 99.90 0 5000 10000 12000 0 0 5000 10000 u", 0, false, false},       // Another primary
            },
        },
        {
            "Keeping only bestN = 2 secondary overlaps. The unit under test does not sort, just picks the first overlaps it finds",
            {
                {"000000000 000000001 -1000 99.90 0 5000 10000 12000 0 0 5000 10000 u", 0, false, false},       // Primary
                {"000000000 000000001 -1000 99.90 0 10000 12000 12000 0 7000 9000 10000 u", 0, true, false},    // Supplementary
                {"000000000 000000002 -1000 99.90 0 5000 10000 10000 0 0 5000 10000 u", 1, false, false},       // Secondary
                {"", 0, false, true},
                {"000000000 000000003 -1000 99.90 0 5000 10000 10000 0 0 5000 10000 u", 1, false, false},       // Secondary
                {"000000000 000000004 -1000 99.90 0 5000 10000 10000 0 0 5000 10000 u", 1, false, false},       // Secondary
                {"000000005 000000006 -1000 99.90 0 5000 10000 12000 0 0 5000 10000 u", 0, false, false},       // Another primary
            },
            2,
            {
                {"000000000 000000001 -1000 99.90 0 5000 10000 12000 0 0 5000 10000 u", 0, false, false},       // Primary
                {"000000000 000000001 -1000 99.90 0 10000 12000 12000 0 7000 9000 10000 u", 0, true, false},    // Supplementary
                {"000000000 000000002 -1000 99.90 0 5000 10000 10000 0 0 5000 10000 u", 1, false, false},       // Secondary
                {"000000000 000000003 -1000 99.90 0 5000 10000 10000 0 0 5000 10000 u", 1, false, false},       // Secondary
                {"000000005 000000006 -1000 99.90 0 5000 10000 12000 0 0 5000 10000 u", 0, false, false},       // Another primary
            },
        },
        {
            "Keeping all secondary overlaps regardless of how many there are (bestN = -1).",
            {
                {"000000000 000000001 -1000 99.90 0 5000 10000 12000 0 0 5000 10000 u", 0, false, false},       // Primary
                {"000000000 000000001 -1000 99.90 0 10000 12000 12000 0 7000 9000 10000 u", 0, true, false},    // Supplementary
                {"000000000 000000002 -1000 99.90 0 5000 10000 10000 0 0 5000 10000 u", 1, false, false},       // Secondary
                {"", 0, false, true},
                {"000000000 000000003 -1000 99.90 0 5000 10000 10000 0 0 5000 10000 u", 1, false, false},       // Secondary
                {"000000000 000000004 -1000 99.90 0 5000 10000 10000 0 0 5000 10000 u", 1, false, false},       // Secondary
                {"000000005 000000006 -1000 99.90 0 5000 10000 12000 0 0 5000 10000 u", 0, false, false},       // Another primary
            },
            -1,
            {
                {"000000000 000000001 -1000 99.90 0 5000 10000 12000 0 0 5000 10000 u", 0, false, false},       // Primary
                {"000000000 000000001 -1000 99.90 0 10000 12000 12000 0 7000 9000 10000 u", 0, true, false},    // Supplementary
                {"000000000 000000002 -1000 99.90 0 5000 10000 10000 0 0 5000 10000 u", 1, false, false},       // Secondary
                {"000000000 000000003 -1000 99.90 0 5000 10000 10000 0 0 5000 10000 u", 1, false, false},       // Secondary
                {"000000000 000000004 -1000 99.90 0 5000 10000 10000 0 0 5000 10000 u", 1, false, false},       // Secondary
                {"000000005 000000006 -1000 99.90 0 5000 10000 12000 0 0 5000 10000 u", 0, false, false},       // Another primary
            },
        },
        {
            "Test that the condensing will work fine when the first overlap is nullptr. Keeping all secondary overlaps regardless of how many there are (bestN = -1).",
            {
                {"", 0, false, true},                                                                      // Nullptr mapping.
                {"000000000 000000003 7630 99.77 0 0 3845 11811 1 0 3838 10150 *", 1, false, false},       // Primary
                {"000000000 000000002 9854 99.90 0 0 4946 11811 0 4181 9125 9139 *", 1, false, false},       // Secondary
                {"000000000 000000004 11982 99.82 0 0 6031 11811 0 2564 8590 8590 *", 1, false, false},       // Secondary
                {"000000000 000000005 15756 99.94 0 0 7895 11811 1 0 7894 11307 *", 1, false, false},       // Secondary
                {"000000000 000000001 18014 99.74 0 0 9077 11811 0 981 10059 10059 *", 1, false, false},       // Secondary
            },
            -1,
            {
                {"000000000 000000003 7630 99.77 0 0 3845 11811 1 0 3838 10150 *", 1, false, false},       // Primary
                {"000000000 000000002 9854 99.90 0 0 4946 11811 0 4181 9125 9139 *", 1, false, false},       // Secondary
                {"000000000 000000004 11982 99.82 0 0 6031 11811 0 2564 8590 8590 *", 1, false, false},       // Secondary
                {"000000000 000000005 15756 99.94 0 0 7895 11811 1 0 7894 11307 *", 1, false, false},       // Secondary
                {"000000000 000000001 18014 99.74 0 0 9077 11811 0 981 10059 10059 *", 1, false, false},       // Secondary
            },
        },
    };
    // clang-format on

    PacBio::Pancake::MapperCLRSettings settings;
    settings.map.freqPercentile = 0.000;

    for (const auto& data : testData) {
        SCOPED_TRACE(data.testName);

        std::cerr << "testName = " << data.testName << "\n";

        // Prepare the input for the test.
        std::vector<std::unique_ptr<ChainedRegion>> mappings;
        for (const auto& vals : data.inMappings) {
            const std::string& ovlStr = std::get<0>(vals);
            const int32_t priority = std::get<1>(vals);
            const bool isSupplementary = std::get<2>(vals);
            const bool isOverlapNullptr = std::get<3>(vals);

            OverlapPtr ovl =
                (isOverlapNullptr) ? nullptr : PacBio::Pancake::ParseM4OverlapFromString(ovlStr);
            auto newChainedRegion = std::make_unique<ChainedRegion>();
            newChainedRegion->mapping = std::move(ovl);
            newChainedRegion->priority = priority;
            newChainedRegion->isSupplementary = isSupplementary;
            mappings.emplace_back(std::move(newChainedRegion));
        }

        // Run the unit under test.
        CondenseMappings(mappings, data.bestNSecondary);

        // Convert the results to a form we can compare.
        std::vector<std::tuple<std::string, int32_t, bool, bool>> results;
        for (const auto& mapping : mappings) {
            if (mapping->mapping == nullptr) {
                results.emplace_back(
                    std::make_tuple(nullptr, mapping->priority, mapping->isSupplementary, true));
            } else {
                results.emplace_back(std::make_tuple(
                    OverlapWriterBase::PrintOverlapAsM4(*mapping->mapping, "", "", true, false),
                    mapping->priority, mapping->isSupplementary, false));
            }
        }

        EXPECT_EQ(data.expected, results);
    }
}

TEST(MapperCLR, CheckSelfHitPolicyAndSkippingSymmetrical)
{
    struct TestData
    {
        std::string testName;
        std::string targetFile;
        std::string queryFile;
        PacBio::Pancake::MapperCLRSettings settings;
        std::string expectedOverlapsPath;
    };

    PacBio::Pancake::MapperCLRSettings settingsDefaultPolicy;
    {
        auto& settings = settingsDefaultPolicy;
        settings.map.bestNSecondary = -1;
        settings.map.secondaryMinScoreFraction = 0.0;
        settings.map.secondaryAllowedOverlapFractionTarget = 0.0;
        settings.map.secondaryAllowedOverlapFractionQuery = 0.0;
        settings.map.skipSymmetricOverlaps = false;
        settings.map.selfHitPolicy = PacBio::Pancake::MapperSelfHitPolicy::DEFAULT;
        settings.align.selfHitPolicy = PacBio::Pancake::MapperSelfHitPolicy::DEFAULT;
        settings.map.freqPercentile = 0.000;
        settings.map.seedParams =
            PacBio::Pancake::SeedDB::SeedDBParameters{19, 10, 0, false, false, 255, true};
        settings.map.seedParamsFallback =
            PacBio::Pancake::SeedDB::SeedDBParameters{10, 5, 0, false, false, 255, true};
    }

    PacBio::Pancake::MapperCLRSettings settingsSkipSelfHitsInBothMapAndAlign;
    {
        auto& settings = settingsSkipSelfHitsInBothMapAndAlign;
        settings.map.bestNSecondary = -1;
        settings.map.secondaryMinScoreFraction = 0.0;
        settings.map.secondaryAllowedOverlapFractionTarget = 0.0;
        settings.map.secondaryAllowedOverlapFractionQuery = 0.0;
        settings.map.skipSymmetricOverlaps = false;
        settings.map.selfHitPolicy = PacBio::Pancake::MapperSelfHitPolicy::SKIP;
        settings.align.selfHitPolicy = PacBio::Pancake::MapperSelfHitPolicy::SKIP;
        settings.map.freqPercentile = 0.000;
        settings.map.seedParams =
            PacBio::Pancake::SeedDB::SeedDBParameters{19, 10, 0, false, false, 255, true};
        settings.map.seedParamsFallback =
            PacBio::Pancake::SeedDB::SeedDBParameters{10, 5, 0, false, false, 255, true};
    }

    PacBio::Pancake::MapperCLRSettings settingsPerfectAlignSelfHitsInBothMapAndAlign;
    {
        auto& settings = settingsPerfectAlignSelfHitsInBothMapAndAlign;
        settings.map.bestNSecondary = -1;
        settings.map.secondaryMinScoreFraction = 0.0;
        settings.map.secondaryAllowedOverlapFractionTarget = 0.0;
        settings.map.secondaryAllowedOverlapFractionQuery = 0.0;
        settings.map.skipSymmetricOverlaps = false;
        settings.map.selfHitPolicy = PacBio::Pancake::MapperSelfHitPolicy::PERFECT_ALIGNMENT;
        settings.align.selfHitPolicy = PacBio::Pancake::MapperSelfHitPolicy::PERFECT_ALIGNMENT;
        settings.map.freqPercentile = 0.000;
        settings.map.seedParams =
            PacBio::Pancake::SeedDB::SeedDBParameters{19, 10, 0, false, false, 255, true};
        settings.map.seedParamsFallback =
            PacBio::Pancake::SeedDB::SeedDBParameters{10, 5, 0, false, false, 255, true};
    }

    PacBio::Pancake::MapperCLRSettings settingsSkipSymmetricOverlaps;
    {
        auto& settings = settingsSkipSymmetricOverlaps;
        settings.map.bestNSecondary = -1;
        settings.map.secondaryMinScoreFraction = 0.0;
        settings.map.secondaryAllowedOverlapFractionTarget = 0.0;
        settings.map.secondaryAllowedOverlapFractionQuery = 0.0;
        settings.map.skipSymmetricOverlaps = true;
        settings.map.selfHitPolicy = PacBio::Pancake::MapperSelfHitPolicy::DEFAULT;
        settings.align.selfHitPolicy = PacBio::Pancake::MapperSelfHitPolicy::DEFAULT;
        settings.map.freqPercentile = 0.000;
        settings.map.seedParams =
            PacBio::Pancake::SeedDB::SeedDBParameters{19, 10, 0, false, false, 255, true};
        settings.map.seedParamsFallback =
            PacBio::Pancake::SeedDB::SeedDBParameters{10, 5, 0, false, false, 255, true};
    }

    PacBio::Pancake::MapperCLRSettings settingsSkipSelfAndSymmetricOverlaps;
    {
        auto& settings = settingsSkipSelfAndSymmetricOverlaps;
        settings.map.bestNSecondary = -1;
        settings.map.secondaryMinScoreFraction = 0.0;
        settings.map.secondaryAllowedOverlapFractionTarget = 0.0;
        settings.map.secondaryAllowedOverlapFractionQuery = 0.0;
        settings.map.skipSymmetricOverlaps = true;
        settings.map.selfHitPolicy = PacBio::Pancake::MapperSelfHitPolicy::SKIP;
        settings.align.selfHitPolicy = PacBio::Pancake::MapperSelfHitPolicy::SKIP;
        settings.map.freqPercentile = 0.000;
        settings.map.seedParams =
            PacBio::Pancake::SeedDB::SeedDBParameters{19, 10, 0, false, false, 255, true};
        settings.map.seedParamsFallback =
            PacBio::Pancake::SeedDB::SeedDBParameters{10, 5, 0, false, false, 255, true};
    }

    PacBio::Pancake::MapperCLRSettings settingsSkipSelfInMappingButDefaultInAlignment;
    {
        auto& settings = settingsSkipSelfInMappingButDefaultInAlignment;
        settings.map.bestNSecondary = -1;
        settings.map.secondaryMinScoreFraction = 0.0;
        settings.map.secondaryAllowedOverlapFractionTarget = 0.0;
        settings.map.secondaryAllowedOverlapFractionQuery = 0.0;
        settings.map.skipSymmetricOverlaps = false;
        settings.map.selfHitPolicy = PacBio::Pancake::MapperSelfHitPolicy::SKIP;
        settings.align.selfHitPolicy = PacBio::Pancake::MapperSelfHitPolicy::DEFAULT;
        settings.map.freqPercentile = 0.000;
        settings.map.seedParams =
            PacBio::Pancake::SeedDB::SeedDBParameters{19, 10, 0, false, false, 255, true};
        settings.map.seedParamsFallback =
            PacBio::Pancake::SeedDB::SeedDBParameters{10, 5, 0, false, false, 255, true};
    }

    PacBio::Pancake::MapperCLRSettings settingsDefaultSelfInMappingButSkipInAlignment;
    {
        auto& settings = settingsDefaultSelfInMappingButSkipInAlignment;
        settings.map.bestNSecondary = -1;
        settings.map.secondaryMinScoreFraction = 0.0;
        settings.map.secondaryAllowedOverlapFractionTarget = 0.0;
        settings.map.secondaryAllowedOverlapFractionQuery = 0.0;
        settings.map.skipSymmetricOverlaps = false;
        settings.map.selfHitPolicy = PacBio::Pancake::MapperSelfHitPolicy::DEFAULT;
        settings.align.selfHitPolicy = PacBio::Pancake::MapperSelfHitPolicy::SKIP;
        settings.map.freqPercentile = 0.000;
        settings.map.seedParams =
            PacBio::Pancake::SeedDB::SeedDBParameters{19, 10, 0, false, false, 255, true};
        settings.map.seedParamsFallback =
            PacBio::Pancake::SeedDB::SeedDBParameters{10, 5, 0, false, false, 255, true};
    }

    PacBio::Pancake::MapperCLRSettings settingsMockSelfInMappingButDefaultInAlignment;
    {
        auto& settings = settingsMockSelfInMappingButDefaultInAlignment;
        settings.map.bestNSecondary = -1;
        settings.map.secondaryMinScoreFraction = 0.0;
        settings.map.secondaryAllowedOverlapFractionTarget = 0.0;
        settings.map.secondaryAllowedOverlapFractionQuery = 0.0;
        settings.map.skipSymmetricOverlaps = false;
        settings.map.selfHitPolicy = PacBio::Pancake::MapperSelfHitPolicy::PERFECT_ALIGNMENT;
        settings.align.selfHitPolicy = PacBio::Pancake::MapperSelfHitPolicy::DEFAULT;
        settings.map.freqPercentile = 0.000;
        settings.map.seedParams =
            PacBio::Pancake::SeedDB::SeedDBParameters{19, 10, 0, false, false, 255, true};
        settings.map.seedParamsFallback =
            PacBio::Pancake::SeedDB::SeedDBParameters{10, 5, 0, false, false, 255, true};
    }

    PacBio::Pancake::MapperCLRSettings settingsDefaultSelfInMappingButMockInAlignment;
    {
        auto& settings = settingsDefaultSelfInMappingButMockInAlignment;
        settings.map.bestNSecondary = -1;
        settings.map.secondaryMinScoreFraction = 0.0;
        settings.map.secondaryAllowedOverlapFractionTarget = 0.0;
        settings.map.secondaryAllowedOverlapFractionQuery = 0.0;
        settings.map.skipSymmetricOverlaps = false;
        settings.map.selfHitPolicy = PacBio::Pancake::MapperSelfHitPolicy::DEFAULT;
        settings.align.selfHitPolicy = PacBio::Pancake::MapperSelfHitPolicy::PERFECT_ALIGNMENT;
        settings.map.freqPercentile = 0.000;
        settings.map.seedParams =
            PacBio::Pancake::SeedDB::SeedDBParameters{19, 10, 0, false, false, 255, true};
        settings.map.seedParamsFallback =
            PacBio::Pancake::SeedDB::SeedDBParameters{10, 5, 0, false, false, 255, true};
    }

    // clang-format off
    std::vector<TestData> testData = {
        {
            "Overlap the same set of reads with itself.",
            // Target.
            PacBio::PancakeTestsConfig::Data_Dir + "/ovl-clr/reads.pile1-5prime.fasta",
            // Query.
            PacBio::PancakeTestsConfig::Data_Dir + "/ovl-clr/reads.pile1-5prime.fasta",
            settingsDefaultPolicy,
            // Expected overlaps.
            PacBio::PancakeTestsConfig::Data_Dir + "/ovl-clr/reads.pile1-5prime.out.all_vs_all.ksw2.m4",
        },
        {
            "Skip self hits.",
            // Target.
            PacBio::PancakeTestsConfig::Data_Dir + "/ovl-clr/reads.pile1-5prime.fasta",
            // Query.
            PacBio::PancakeTestsConfig::Data_Dir + "/ovl-clr/reads.pile1-5prime.fasta",
            settingsSkipSelfHitsInBothMapAndAlign,
            // Expected overlaps.
            PacBio::PancakeTestsConfig::Data_Dir + "/ovl-clr/reads.pile1-5prime.out.all_vs_all_no_self_hits.ksw2.m4",
        },
        {
            "Mock perfect overlaps",
            // Target.
            PacBio::PancakeTestsConfig::Data_Dir + "/ovl-clr/reads.pile1-5prime.fasta",
            // Query.
            PacBio::PancakeTestsConfig::Data_Dir + "/ovl-clr/reads.pile1-5prime.fasta",
            settingsPerfectAlignSelfHitsInBothMapAndAlign,
            // Expected overlaps.
            PacBio::PancakeTestsConfig::Data_Dir + "/ovl-clr/reads.pile1-5prime.out.all_vs_all.ksw2.m4",
        },
        {
            "Skip symmetric overlaps.",
            // Target.
            PacBio::PancakeTestsConfig::Data_Dir + "/ovl-clr/reads.pile1-5prime.fasta",
            // Query.
            PacBio::PancakeTestsConfig::Data_Dir + "/ovl-clr/reads.pile1-5prime.fasta",
            settingsSkipSymmetricOverlaps,
            // Expected overlaps.
            PacBio::PancakeTestsConfig::Data_Dir + "/ovl-clr/reads.pile1-5prime.out.all_vs_all_no_symmetric.ksw2.m4",
        },
        {
            "Skip self and symmetric overlaps.",
            // Target.
            PacBio::PancakeTestsConfig::Data_Dir + "/ovl-clr/reads.pile1-5prime.fasta",
            // Query.
            PacBio::PancakeTestsConfig::Data_Dir + "/ovl-clr/reads.pile1-5prime.fasta",
            settingsSkipSelfAndSymmetricOverlaps,
            // Expected overlaps.
            PacBio::PancakeTestsConfig::Data_Dir + "/ovl-clr/reads.pile1-5prime.out.all_vs_all_no_self_hits_no_symmetric.ksw2.m4",
        },
        {
            "Skip self hits in the mapping stage, but use the default policy during alignment. This should skip the self hits completely.",
            // Target.
            PacBio::PancakeTestsConfig::Data_Dir + "/ovl-clr/reads.pile1-5prime.fasta",
            // Query.
            PacBio::PancakeTestsConfig::Data_Dir + "/ovl-clr/reads.pile1-5prime.fasta",
            settingsSkipSelfInMappingButDefaultInAlignment,
            // Expected overlaps.
            PacBio::PancakeTestsConfig::Data_Dir + "/ovl-clr/reads.pile1-5prime.out.all_vs_all_no_self_hits.ksw2.m4",
        },
        {
            "Skip self hits in the alignment stage, but use the default policy during mapping. This should skip the self hits completely.",
            // Target.
            PacBio::PancakeTestsConfig::Data_Dir + "/ovl-clr/reads.pile1-5prime.fasta",
            // Query.
            PacBio::PancakeTestsConfig::Data_Dir + "/ovl-clr/reads.pile1-5prime.fasta",
            settingsDefaultSelfInMappingButSkipInAlignment,
            // Expected overlaps.
            PacBio::PancakeTestsConfig::Data_Dir + "/ovl-clr/reads.pile1-5prime.out.all_vs_all_no_self_hits.ksw2.m4",
        },
        {
            "Mock perfect overlaps in the mapping stage, but use the default policy during alignment. This should report proper alignments, like everything was default.",
            // Target.
            PacBio::PancakeTestsConfig::Data_Dir + "/ovl-clr/reads.pile1-5prime.fasta",
            // Query.
            PacBio::PancakeTestsConfig::Data_Dir + "/ovl-clr/reads.pile1-5prime.fasta",
            settingsMockSelfInMappingButDefaultInAlignment,
            // Expected overlaps.
            PacBio::PancakeTestsConfig::Data_Dir + "/ovl-clr/reads.pile1-5prime.out.all_vs_all.ksw2.m4",
        },
        {
            "Mock perfect overlaps in the alignment stage, but use the default policy during mapping. This should report proper alignments, like everything was default.",
            // Target.
            PacBio::PancakeTestsConfig::Data_Dir + "/ovl-clr/reads.pile1-5prime.fasta",
            // Query.
            PacBio::PancakeTestsConfig::Data_Dir + "/ovl-clr/reads.pile1-5prime.fasta",
            settingsDefaultSelfInMappingButMockInAlignment,
            // Expected overlaps.
            PacBio::PancakeTestsConfig::Data_Dir + "/ovl-clr/reads.pile1-5prime.out.all_vs_all.ksw2.m4",
        },
    };
    // clang-format on

    auto CreateFastaSequenceCached = [](const std::vector<PacBio::BAM::FastaSequence>& seqs)
        -> std::vector<PacBio::Pancake::FastaSequenceCached> {
        std::vector<PacBio::Pancake::FastaSequenceCached> records;
        for (size_t id = 0; id < seqs.size(); ++id) {
            const auto& record = seqs[id];
            records.emplace_back(PacBio::Pancake::FastaSequenceCached(
                record.Name(), record.Bases().c_str(), record.Bases().size(), id));
        }
        return records;
    };

    for (const auto& data : testData) {
        // Load the sequences from files, and take only the first one.
        const std::vector<PacBio::BAM::FastaSequence> allTargetSeqs =
            PacBio::PancakeTests::HelperLoadFasta(data.targetFile);
        const std::vector<PacBio::BAM::FastaSequence> allQuerySeqs =
            PacBio::PancakeTests::HelperLoadFasta(data.queryFile);
        const std::vector<std::string> expectedOverlaps =
            PacBio::PancakeTests::HelperLoadFile(data.expectedOverlapsPath);

        const std::vector<PacBio::Pancake::FastaSequenceCached> targets =
            CreateFastaSequenceCached(allTargetSeqs);
        const std::vector<PacBio::Pancake::FastaSequenceCached> queries =
            CreateFastaSequenceCached(allQuerySeqs);

        SCOPED_TRACE(data.testName);

        std::cerr << "testName = " << data.testName << "\n";

        // Run mapping.
        PacBio::Pancake::MapperCLR mapper(data.settings);
        std::vector<PacBio::Pancake::MapperBaseResult> result =
            mapper.MapAndAlign(targets, queries);

        // Convert results to strings so we can compare easily.
        std::vector<std::string> resultsStr;
        for (const auto& queryMappings : result) {
            for (const auto& mapping : queryMappings.mappings) {
                resultsStr.emplace_back(PacBio::Pancake::OverlapWriterBase::PrintOverlapAsM4(
                    *mapping->mapping, "", "", true, false));
            }
        }

        // Sort the output so that we can compare sets easily.
        auto expected = expectedOverlaps;
        std::sort(expected.begin(), expected.end());
        std::sort(resultsStr.begin(), resultsStr.end());

        // std::cerr << "Expected:\n";
        // for (const auto& ovlStr : expected) {
        //     std::cerr << ovlStr << "\n";
        // }
        // std::cerr << "Results:\n";
        // for (const auto& ovlStr : resultsStr) {
        //     std::cerr << ovlStr << "\n";
        // }

        EXPECT_EQ(expected, resultsStr);
    }
    // exit(1);
}
