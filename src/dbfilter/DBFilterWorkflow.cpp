// Authors: Ivan Sovic

#include "DBFilterWorkflow.h"

#include <algorithm>
#include <iostream>
#include <random>
#include <string>
#include <unordered_set>
#include <vector>

#include <pbcopper/logging/LogLevel.h>
#include <pbcopper/logging/Logging.h>

#include <boost/algorithm/string/predicate.hpp>

#include <pacbio/seeddb/SeedDBIndexCache.h>
#include <pacbio/seqdb/SeqDBIndexCache.h>
#include <pacbio/seqdb/SeqDBReader.h>
#include <pacbio/seqdb/SeqDBWriter.h>
#include <pacbio/seqdb/Util.h>
#include <pacbio/util/TicToc.h>

#include "DBFilterSettings.h"

namespace PacBio {
namespace Pancake {

void PerformSeqDBSequenceLineSampling(std::vector<SeqDBSequenceLine>& outSeqLines,
                                      const std::vector<SeqDBSequenceLine>& inSeqLines,
                                      const SamplingType& sampling, int64_t sampledBases,
                                      const int64_t randomSeed,
                                      const std::unordered_set<std::string>& filterList,
                                      const FilterListType& filterType)
{
    outSeqLines.size();

    auto CheckFilterShouldKeep = [&](const std::string& header) {
        if (filterType == FilterListType::Blacklist &&
            filterList.find(header) != filterList.end()) {
            return false;
        }
        if (filterType == FilterListType::Whitelist &&
            filterList.find(header) == filterList.end()) {
            return false;
        }
        return true;
    };

    if (sampling == SamplingType::Linear) {
        int64_t totalBases = 0;
        for (int32_t lastLine = 0;
             lastLine < static_cast<int32_t>(inSeqLines.size()) && totalBases < sampledBases;
             ++lastLine) {
            const auto& sl = inSeqLines[lastLine];
            // Filter sequences.
            if (CheckFilterShouldKeep(sl.header) == false) {
                continue;
            }
            totalBases += sl.numBases;
            outSeqLines.emplace_back(sl);
            if (totalBases >= sampledBases) {
                break;
            }
        }

    } else if (sampling == SamplingType::Random) {
        std::random_device rd;
        const uint64_t seed =
            (randomSeed < 0) ? std::mt19937::default_seed : static_cast<uint64_t>(randomSeed);
        std::mt19937 eng(seed);
        if (randomSeed < 0) {
            eng = std::mt19937(rd());
        }

        // Shuffle the permutation.
        std::vector<int32_t> permutation(inSeqLines.size());
        std::iota(permutation.begin(), permutation.end(), 0);
        for (size_t i = 0; i < permutation.size(); ++i) {
            size_t j = eng() % permutation.size();
            std::swap(permutation[i], permutation[j]);
        }

        // The rest is similar to Linear sampling, but with an index redirection.
        int64_t totalBases = 0;
        for (size_t i = 0; i < permutation.size() && totalBases < sampledBases; ++i) {
            const auto& sl = inSeqLines[permutation[i]];
            // Filter sequences.
            if (CheckFilterShouldKeep(sl.header) == false) {
                continue;
            }
            totalBases += sl.numBases;
            outSeqLines.emplace_back(sl);
            if (totalBases >= sampledBases) {
                break;
            }
        }
        // Sort by sequence ID, to preserve the cache coherency if possible.
        std::sort(outSeqLines.begin(), outSeqLines.end(),
                  [](const auto& a, const auto& b) { return a.seqId < b.seqId; });

    } else if (sampling == SamplingType::None) {
        for (size_t i = 0; i < inSeqLines.size(); ++i) {
            const auto& sl = inSeqLines[i];
            // Filter sequences.
            if (CheckFilterShouldKeep(sl.header) == false) {
                continue;
            }
            outSeqLines.emplace_back(sl);
        }

    } else {
        throw std::runtime_error("Unknown sampling method!");
    }
}

std::unordered_set<std::string> ParseFilterList(const std::string& filterListPath)
{
    std::unordered_set<std::string> filterList;
    std::ifstream ifs(filterListPath);
    if (ifs.is_open() == false) {
        throw std::runtime_error("Could not open blacklist file '" + filterListPath + "'!");
    }
    std::string line;
    while (std::getline(ifs, line)) {
        filterList.emplace(line);
    }
    return filterList;
}

int DBFilterWorkflow::Runner(const PacBio::CLI_v2::Results& options)
{
    DBFilterSettings settings{options};

    std::string inSeqDBFile = settings.InputPrefix + ".seqdb";
    // std::string inSeedDBFile = settings.InputPrefix + ".seeddb";
    std::string outSeqDBFile = settings.OutputPrefix + ".seqdb";
    // std::string outSeedDBFile = settings.OutputPrefix + ".seeddb";

    PBLOG_INFO << "Loading the input DBs.";
    TicToc ttLoad;

    // Load the target SeqDB caches.
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> inSeqDBCache =
        PacBio::Pancake::LoadSeqDBIndexCache(inSeqDBFile);
    PBLOG_INFO << "After loading target seq cache: " << ttLoad.VerboseSecs(true);

    // // Load the target SeedDB caches.
    // std::shared_ptr<PacBio::Pancake::SeedDBIndexCache> inSeedDBCache =
    //     PacBio::Pancake::LoadSeedDBIndexCache(inSeedDBFile);
    // PBLOG_INFO << "After loading target seed cache: " << ttLoad.VerboseSecs(true);
    ttLoad.Stop();

    std::unordered_set<std::string> blacklist;
    if (settings.FilterListPath.empty() == false) {
        blacklist = ParseFilterList(settings.FilterListPath);
    }
    PBLOG_INFO << "Filter list size: " << blacklist.size();

    // Construct a filtered SeqDB cache.
    // Do not normalize yet, we might need the original sequence IDs for fetching.
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> filteredSeqDBCache =
        std::make_shared<PacBio::Pancake::SeqDBIndexCache>();
    filteredSeqDBCache->indexFilename = outSeqDBFile;
    SplitPath(filteredSeqDBCache->indexFilename, filteredSeqDBCache->indexParentFolder,
              filteredSeqDBCache->indexBasename);
    filteredSeqDBCache->version = inSeqDBCache->version;
    filteredSeqDBCache->compressionLevel = inSeqDBCache->compressionLevel;
    // The data will not be copied (only the index), so the file lines are the same.
    filteredSeqDBCache->fileLines = inSeqDBCache->fileLines;
    PerformSeqDBSequenceLineSampling(filteredSeqDBCache->seqLines, inSeqDBCache->seqLines,
                                     settings.Sampling, settings.SampleBases, settings.RandomSeed,
                                     blacklist, settings.FilterType);

    PBLOG_INFO << "Original SeqDB sequences: " << inSeqDBCache->seqLines.size();
    PBLOG_INFO << "Filtered SeqDB sequences: " << filteredSeqDBCache->seqLines.size();

    // If the Consolidate option is not specified, then we only need to modify the
    // index cache. Otherwise, we need to create a copy of the data using the Writer.
    if (settings.Consolidate) {
        // Create a reader so we can fetch sequences.
        PacBio::Pancake::SeqDBReader reader(inSeqDBCache);

        // Create a new writer, which uses the actual output prefix.
        auto writer = PacBio::Pancake::CreateSeqDBWriter(
            settings.OutputPrefix, settings.CompressionLevel, settings.BufferSize,
            settings.BlockSize, settings.SplitBlocks);

        // Fetch and copy all the sequences into the new DB.
        Pancake::FastaSequenceId record;
        for (const auto& sl : filteredSeqDBCache->seqLines) {
            reader.GetSequence(record, sl.seqId);
            writer->AddSequence(record.Name(), record.Bases());
        }

    } else {
        NormalizeSeqDBIndexCache(*filteredSeqDBCache, settings.BlockSize);
        std::unique_ptr<FILE, FileDeleter> fpOutSeqDBCache =
            PacBio::Pancake::OpenFile(outSeqDBFile.c_str(), "w");
        WriteSeqDBIndexCache(fpOutSeqDBCache.get(), *filteredSeqDBCache);
    }

    PBLOG_INFO << "Done writing of the filtered DB.";

    return EXIT_SUCCESS;
}

}  // namespace Pancake
}  // namespace PacBio
