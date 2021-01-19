#include <PancakeTestData.h>
#include <gtest/gtest.h>
#include <pacbio/pancake/MapperBatchCPU.h>
#include <pacbio/pancake/OverlapWriterBase.h>
#include <iostream>
#include "TestHelperUtils.h"

void HelperLoadBatchData(
    const std::vector<std::pair<std::string, std::string>>& batchDataSequenceFiles,
    std::vector<PacBio::Pancake::MapperBatchChunk>& retBatchData,
    std::vector<PacBio::BAM::FastaSequence>& retAllSeqs)
{
    /*
     * This function takes a vector of pairs, where each pair is a pair of filename paths, one
     * for the target sequences and one for the query sequences; and then loads the sequences.
     * One pair is what a single MapperCLR::MapAndAlign would be run on.
     * Here, the MapperBatchCPU will take a vector of those chunks and align them at once.
     *
     * This function returns a flat vector of all sequences that were loaded (targets and queries),
     * and a vector of MapperBatchChunk.
     * The flat vector is needed because MapperBatchChunk uses FastaSequenceCached which uses pointers
     * to those sequences, so the life span of the original data needs to be ensured.
    */

    retBatchData.clear();
    retAllSeqs.clear();

    for (const auto& vals : batchDataSequenceFiles) {
        const auto& targetFile = vals.first;
        const auto& queryFile = vals.second;
        PacBio::Pancake::MapperBatchChunk bd;

        // Load target sequences and construct the FastaSequenceCached objects.
        std::vector<PacBio::BAM::FastaSequence> targetSeqs =
            PacBio::PancakeTests::HelperLoadFasta(targetFile);
        for (size_t seqId = 0; seqId < targetSeqs.size(); ++seqId) {
            const auto& seq = targetSeqs[seqId];
            // const int32_t seqId = retAllSeqs.size();
            retAllSeqs.emplace_back(std::move(seq));
            auto newFsc = PacBio::Pancake::FastaSequenceCached(
                std::to_string(seqId), retAllSeqs.back().Bases().c_str(),
                retAllSeqs.back().Bases().size(), seqId);
            bd.targetSeqs.emplace_back(std::move(newFsc));
        }

        // Load query sequences and construct the FastaSequenceCached objects.
        std::vector<PacBio::BAM::FastaSequence> querySeqs =
            PacBio::PancakeTests::HelperLoadFasta(queryFile);
        for (size_t seqId = 0; seqId < querySeqs.size(); ++seqId) {
            const auto& seq = querySeqs[seqId];
            // const int32_t seqId = retAllSeqs.size();
            retAllSeqs.emplace_back(std::move(seq));
            auto newFsc = PacBio::Pancake::FastaSequenceCached(
                std::to_string(seqId), retAllSeqs.back().Bases().c_str(),
                retAllSeqs.back().Bases().size(), seqId);
            bd.querySeqs.emplace_back(std::move(newFsc));
        }

        retBatchData.emplace_back(std::move(bd));
    }
}

std::vector<std::vector<std::string>> HelperFormatBatchMappingResults(
    const std::vector<std::vector<PacBio::Pancake::MapperBaseResult>>& results)
{
    std::vector<std::vector<std::string>> resultsStr;
    for (size_t i = 0; i < results.size(); ++i) {
        const auto& resultsForBatchElement = results[i];
        std::vector<std::string> chunkResults;
        for (size_t j = 0; j < resultsForBatchElement.size(); ++j) {
            const auto& queryMappings = resultsForBatchElement[j];
            for (const auto& mapping : queryMappings.mappings) {
                // std::cerr << PacBio::Pancake::OverlapWriterBase::PrintOverlapAsM4(
                //                  mapping->mapping, "", "", true, false)
                //           << "\n";

                chunkResults.emplace_back(PacBio::Pancake::OverlapWriterBase::PrintOverlapAsM4(
                    mapping->mapping, "", "", true, false));
            }
        }
        resultsStr.emplace_back(std::move(chunkResults));
    }
    return resultsStr;
}

TEST(MapperBatchCPU, BatchMapping_ArrayOfTests)
{
    using namespace PacBio::Pancake;

    struct TestData
    {
        std::string testName;
        std::vector<std::pair<std::string, std::string>> batchData;
        PacBio::Pancake::SeedDB::SeedDBParameters seedParamsPrimary;
        PacBio::Pancake::SeedDB::SeedDBParameters seedParamsFallback;
        std::vector<std::vector<std::string>> expectedOverlaps;
    };

    // clang-format off
    std::vector<TestData> testData = {
        {
            "Batch of multiple query/target vectors.",
            {
                {
                    PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-2-real-insert-target.fasta",
                    PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-2-real-insert-query.fasta",
                },
                {
                    PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-3-real-subseq-full-global-match-target.fasta",
                    PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-3-real-subseq-full-global-match-query.fasta",
                },
                {
                    PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-4-small-revcmp-perfect-aln-target.fasta",
                    PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-4-small-revcmp-perfect-aln-query.fasta",
                },
                {
                    PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-5-revcmp-aln-target.fasta",
                    PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-5-revcmp-aln-query.fasta",
                },
                {
                    PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-6-small-revcmp-imperfect-aln-target.fasta",
                    PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-6-small-revcmp-imperfect-aln-query.fasta",
                },
                {
                    PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-1-poor-aln-overlapping-seeds.target.fasta",
                    PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-1-poor-aln-overlapping-seeds.query.fasta",
                },

            },
            // SeedParams - primary.
            PacBio::Pancake::SeedDB::SeedDBParameters{19, 10, 0, false, false, 255, true},
            // SeedParams - fallback.
            PacBio::Pancake::SeedDB::SeedDBParameters{10, 5, 0, false, false, 255, true},
            // Expected results.
            {
                {
                    "000000000 000000000 -16391 81.12 0 0 18779 18779 0 0 18864 18865 *",
                },
                {
                    "000000000 000000000 -8049 98.71 0 0 8111 8111 0 0 8138 8138 *"
                },
                {
                    "000000000 000000000 -150 100.00 0 0 150 150 1 0 150 150 *"
                },
                {
                    "000000000 000000000 -17688 74.58 0 0 21382 22015 1 701 22028 22028 *"
                },
                {
                    "000000000 000000000 -150 96.77 0 0 155 155 1 0 150 150 *"
                },
                {
                    "000000000 000000000 -1345 68.58 0 6209 7938 43446 0 7261 8999 46238 *"
                },

            },
        },
    };
    // clang-format on

    PacBio::Pancake::MapperCLRSettings settings;
    settings.freqPercentile = 0.000;

    for (const auto& data : testData) {
        // Debug info.
        SCOPED_TRACE(data.testName);
        std::cerr << "testName = " << data.testName << "\n";

        // Load the batch sequence data. The helper function takes
        // a vector of target-query filename pairs.
        std::vector<MapperBatchChunk> batchData;
        std::vector<PacBio::BAM::FastaSequence> allSeqs;
        HelperLoadBatchData(data.batchData, batchData, allSeqs);

        // Set the seed parameter settings and create a mapper.
        settings.seedParams = data.seedParamsPrimary;
        settings.seedParamsFallback = data.seedParamsFallback;
        PacBio::Pancake::MapperBatchCPU mapper(settings, 1);

        // Run the unit under test.
        // std::vector<std::vector<MapperBaseResult>> results = mapper.DummyMapAndAlign(batchData);
        std::vector<std::vector<MapperBaseResult>> results = mapper.MapAndAlign(batchData);

        // Format the results for comparison.
        std::vector<std::vector<std::string>> resultsStr = HelperFormatBatchMappingResults(results);

        // Evaluate.
        ASSERT_EQ(data.expectedOverlaps, resultsStr);
    }
}
