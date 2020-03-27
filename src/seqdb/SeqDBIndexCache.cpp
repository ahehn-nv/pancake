// Authors: Ivan Sovic

#include <pacbio/seqdb/SeqDBIndexCache.h>
#include <pacbio/seqdb/Util.h>
#include <array>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>

namespace PacBio {
namespace Pancake {

std::unique_ptr<PacBio::Pancake::SeqDBIndexCache> LoadSeqDBIndexCache(
    const std::string& indexFilename)
{
    FILE* fpIn = fopen(indexFilename.c_str(), "r");
    if (fpIn == NULL) {
        std::ostringstream oss;
        oss << "Could not open file '" << indexFilename << "' for reading!";
        throw std::runtime_error(oss.str());
    }
    auto result = LoadSeqDBIndexCache(fpIn, indexFilename);
    fclose(fpIn);
    return result;
}

std::unique_ptr<PacBio::Pancake::SeqDBIndexCache> LoadSeqDBIndexCache(
    FILE* fpIn, const std::string& indexFilename)
{
    auto cache = std::make_unique<PacBio::Pancake::SeqDBIndexCache>();

    cache->indexFilename = indexFilename;
    SplitPath(indexFilename, cache->indexParentFolder, cache->indexBasename);

    char* line = NULL;
    size_t lineLen = 0;
    ssize_t numRead = 0;
    char buff[2000];  // Maximum string length (file names, headers).

    // Helper function to do cleanup of the C-based allocs.
    auto Cleanup = [&]() {
        if (line) {
            free(line);
            line = NULL;
            lineLen = 0;
        }
    };

    SeqDBFileLine fl;
    SeqDBSequenceLine sl;
    SeqDBBlockLine bl;
    int32_t numRanges = 0;
    int32_t readOffset = 0;
    int32_t offset = 0;
    int32_t numReadItems = 0;
    int32_t totalNumSeqs = 0;

    while ((numRead = getline(&line, &lineLen, fpIn)) != -1) {
        if (lineLen <= 0) {
            Cleanup();
            continue;
        }
        const char token = line[0];
        switch (token) {
            case 'V':
                numReadItems = sscanf(&line[1], "%s", buff);
                cache->version = buff;
                if (numReadItems != 1)
                    throw std::runtime_error("Problem parsing line: '" + std::string(line) + "'.");
                break;
            case 'C':
                numReadItems = sscanf(&line[1], "%d", &(cache->compressionLevel));
                if (numReadItems != 1)
                    throw std::runtime_error("Problem parsing line: '" + std::string(line) + "'.");
                break;
            case 'F':
                numReadItems = sscanf(&line[1], "%d %s %d %lld %lld", &(fl.fileId), buff,
                                      &(fl.numSequences), &(fl.numBytes), &(fl.numCompressedBases));
                if (numReadItems != 5)
                    throw std::runtime_error("Problem parsing line: '" + std::string(line) + "'.");
                fl.filename = buff;
                cache->fileLines.emplace_back(fl);
                totalNumSeqs += fl.numSequences;
                cache->seqLines.reserve(totalNumSeqs);
                break;
            case 'S':
                numReadItems = sscanf(&line[1], "%d %s %d %lld %d %d %d%n", &(sl.seqId), buff,
                                      &(sl.fileId), &(sl.fileOffset), &(sl.numBytes),
                                      &(sl.numBases), &(numRanges), &readOffset);
                if (numReadItems != 7)
                    throw std::runtime_error("Problem parsing line: '" + std::string(line) + "'.");
                if (sl.seqId != static_cast<int32_t>(cache->seqLines.size())) {
                    std::ostringstream oss;
                    oss << "Invalid seqId for line: '" << line
                        << "'. The actual ordinal ID of the sequence line is "
                        << cache->seqLines.size();
                    throw std::runtime_error(oss.str());
                }
                sl.header = buff;
                sl.ranges.clear();
                offset = readOffset + 1;
                for (int32_t i = 0; i < numRanges; ++i) {
                    Range r;
                    numReadItems = sscanf(&line[offset], "%d %d%n", &r.start, &r.end, &readOffset);
                    if (numReadItems != 2)
                        throw std::runtime_error("Problem parsing line: '" + std::string(line) +
                                                 "'.");
                    offset += readOffset + 1;
                    sl.ranges.emplace_back(r);
                }
                cache->seqLines.emplace_back(sl);
                break;
            case 'B':
                numReadItems =
                    sscanf(&line[1], "%d %d %d %lld %lld", &(bl.blockId), &(bl.startSeqId),
                           &(bl.endSeqId), &(bl.numBytes), &(bl.numBases));
                if (numReadItems != 5)
                    throw std::runtime_error("Problem parsing line: '" + std::string(line) + "'.");
                cache->blockLines.emplace_back(bl);
                break;
            default:
                std::ostringstream oss;
                oss << "Unknown token found when parsing the index: " << line[0];
                throw std::runtime_error(oss.str());
                break;
        }
        Cleanup();
    }
    Cleanup();
    return cache;
}

std::unique_ptr<PacBio::Pancake::SeqDBIndexCache> LoadSeqDBIndexCache(
    std::istream& is, const std::string& indexFilename)
{
    auto cache = std::make_unique<PacBio::Pancake::SeqDBIndexCache>();

    cache->indexFilename = indexFilename;
    SplitPath(indexFilename, cache->indexParentFolder, cache->indexBasename);

    std::string line;
    char token;
    while (std::getline(is, line)) {
        if (line.empty()) continue;

        std::istringstream iss(line);
        iss >> token;

        SeqDBFileLine fl;
        SeqDBSequenceLine sl;
        SeqDBBlockLine bl;
        int32_t numRanges = 0;

        switch (token) {
            case 'V':
                iss >> cache->version;
                break;
            case 'C':
                iss >> cache->compressionLevel;
                if (cache->compressionLevel < 0)
                    throw std::runtime_error("Unsupported compression level: " +
                                             std::to_string(cache->compressionLevel));
                break;
            case 'F':
                iss >> fl.fileId >> fl.filename >> fl.numSequences >> fl.numBytes >>
                    fl.numCompressedBases;
                cache->fileLines.emplace_back(fl);
                break;
            case 'S':
                iss >> sl.seqId >> sl.header >> sl.fileId >> sl.fileOffset >> sl.numBytes >>
                    sl.numBases >> numRanges;
                for (int32_t i = 0; i < numRanges; ++i) {
                    Range r;
                    iss >> r.start >> r.end;
                    sl.ranges.emplace_back(r);
                }
                if (sl.seqId != static_cast<int32_t>(cache->seqLines.size())) {
                    std::ostringstream oss;
                    oss << "Invalid seqId for line: '" << line
                        << "'. The actual ordinal ID of the sequence line is "
                        << cache->seqLines.size();
                    throw std::runtime_error(oss.str());
                }
                cache->seqLines.emplace_back(sl);
                break;
            case 'B':
                iss >> bl.blockId >> bl.startSeqId >> bl.endSeqId >> bl.numBytes >> bl.numBases;
                cache->blockLines.emplace_back(bl);
                break;
            default:
                std::ostringstream oss;
                oss << "Unknown token found when parsing the index: " << token;
                throw std::runtime_error(oss.str());
                break;
        }
    }

    if (cache->seqLines.empty())
        throw std::runtime_error("There are no sequences in the input index file: " +
                                 indexFilename);

    return cache;
}

void WriteSeqDBIndexCache(FILE* fpOut, const SeqDBIndexCache& cache)
{
    // An output index file should be open at all times, starting from construction.
    if (fpOut == nullptr) {
        throw std::runtime_error("Cannot write the index because an output file is not open.");
    }

    // Write the version and compression information.
    fprintf(fpOut, "V\t%s\n", cache.version.c_str());
    fprintf(fpOut, "C\t%d\n",
            static_cast<int32_t>(cache.compressionLevel));  // Compression is turned on.

    // Write all the files and their sizes.
    for (const auto& f : cache.fileLines) {
        fprintf(fpOut, "F\t%d\t%s\t%d\t%lld\t%lld\n", f.fileId, f.filename.c_str(), f.numSequences,
                f.numBytes, f.numCompressedBases);
    }

    // Write the indexes of all sequences.
    for (size_t i = 0; i < cache.seqLines.size(); ++i) {
        fprintf(fpOut, "S\t%d\t%s\t%d\t%lld\t%d\t%d", cache.seqLines[i].seqId,
                cache.seqLines[i].header.c_str(), cache.seqLines[i].fileId,
                cache.seqLines[i].fileOffset, cache.seqLines[i].numBytes,
                cache.seqLines[i].numBases);
        fprintf(fpOut, "\t%lu", cache.seqLines[i].ranges.size());
        for (const auto& r : cache.seqLines[i].ranges) {
            fprintf(fpOut, "\t%d\t%d", r.start, r.end);
        }
        fprintf(fpOut, "\n");
    }

    // Write the blocks of all sequences.
    for (size_t i = 0; i < cache.blockLines.size(); ++i) {
        fprintf(fpOut, "B\t%d\t%d\t%d\t%lld\t%lld\n", cache.blockLines[i].blockId,
                cache.blockLines[i].startSeqId, cache.blockLines[i].endSeqId,
                cache.blockLines[i].numBytes, cache.blockLines[i].numBases);
    }
}

void ComputeSeqDBIndexHeaderLookup(const PacBio::Pancake::SeqDBIndexCache& dbCache,
                                   HeaderLookupType& headerToOrdinalId)
{
    headerToOrdinalId.clear();
    headerToOrdinalId.reserve(dbCache.seqLines.size());
    int32_t numRecords = dbCache.seqLines.size();
    for (int32_t i = 0; i < numRecords; ++i) {
        const auto& sl = dbCache.seqLines[i];
        headerToOrdinalId[sl.header] = i;
    }
}

std::vector<SeqDBBlockLine> CreateSeqDBBlocks(const std::vector<SeqDBSequenceLine>& seqLines,
                                              int64_t blockSize)
{
    std::vector<SeqDBBlockLine> blocks;
    int32_t numSeqLines = static_cast<int32_t>(seqLines.size());
    int32_t startId = 0;
    int64_t numBases = 0;
    int64_t numBytes = 0;

    for (int32_t i = 0; i < numSeqLines; ++i) {
        const auto& sl = seqLines[i];

        // Sanity check  that the DB index is valid.
        if (sl.seqId != i) {
            std::ostringstream oss;
            oss << "Invalid SeqDB: sequence ID for '" << sl.header
                << "' is not in line with it's order of appearance in the DB. seqId = " << sl.seqId
                << ", i = " << i;
            throw std::runtime_error(oss.str());
        }

        numBases += sl.numBases;
        numBytes += sl.numBytes;

        // Create a new block when the time is right.
        if (numBases >= blockSize) {
            SeqDBBlockLine bl;
            bl.blockId = blocks.size();
            bl.startSeqId = startId;
            bl.endSeqId = sl.seqId + 1;
            bl.numBytes = numBytes;
            bl.numBases = numBases;
            blocks.emplace_back(bl);
            numBases = 0;
            numBytes = 0;
            startId = sl.seqId + 1;
        }
    }

    // Last block.
    if (startId < numSeqLines) {
        SeqDBBlockLine bl;
        bl.blockId = blocks.size();
        bl.startSeqId = startId;
        bl.endSeqId = numSeqLines;
        bl.numBytes = numBytes;
        bl.numBases = numBases;
        blocks.emplace_back(bl);
    }

    return blocks;
}

void NormalizeSeqDBIndexCache(SeqDBIndexCache& cache, int64_t blockSize)
{
    for (int32_t i = 0; i < static_cast<int32_t>(cache.seqLines.size()); ++i) {
        auto& sl = cache.seqLines[i];
        sl.seqId = i;
    }
    cache.blockLines = CreateSeqDBBlocks(cache.seqLines, blockSize);
}

void SeqDBIndexCache::Validate() const
{
    if (fileLines.empty())
        throw std::runtime_error("There are no file specifications in the input index file.");
    if (seqLines.empty())
        throw std::runtime_error("There are no sequences in the input index file.");
    if (blockLines.empty())
        throw std::runtime_error("There are no blocks in the input index file.");
}

void ValidateSeqDBIndexCache(std::shared_ptr<PacBio::Pancake::SeqDBIndexCache>& indexCache)
{
    // Sanity checks.
    if (indexCache == nullptr) throw std::runtime_error("Provided seqDBCache == nullptr!");
    indexCache->Validate();
}

const SeqDBSequenceLine& SeqDBIndexCache::GetSeqLine(int32_t seqId) const
{
    // Sanity check for the sequence ID.
    if (seqId < 0 || seqId >= static_cast<int32_t>(seqLines.size())) {
        std::ostringstream oss;
        oss << "Invalid seqId. seqId = " << seqId << ", seedLines.size() = " << seqLines.size();
        throw std::runtime_error(oss.str());
    }
    return seqLines[seqId];
}

const SeqDBBlockLine& SeqDBIndexCache::GetBlockLine(int32_t blockId) const
{
    // Sanity check for the sequence ID.
    if (blockId < 0 || blockId >= static_cast<int32_t>(blockLines.size())) {
        std::ostringstream oss;
        oss << "Invalid blockId. blockId = " << blockId
            << ", blockLines.size() = " << blockLines.size();
        throw std::runtime_error(oss.str());
    }
    return blockLines[blockId];
}

const SeqDBFileLine& SeqDBIndexCache::GetFileLine(int32_t fileId) const
{
    // Sanity check.
    if (fileId < 0 || fileId >= static_cast<int32_t>(fileLines.size())) {
        std::ostringstream oss;
        oss << "Invalid fileId. fileId = " << fileId << ", fileLines.size() = " << fileLines.size();
        throw std::runtime_error(oss.str());
    }
    return fileLines[fileId];
}

std::vector<ContiguousFilePart> GetSeqDBContiguousParts(
    const std::shared_ptr<PacBio::Pancake::SeqDBIndexCache>& seqDBIndexCache, int32_t blockId)
{
    const auto& block = seqDBIndexCache->GetBlockLine(blockId);

    if (block.startSeqId < 0 || block.endSeqId < 0 || block.endSeqId < block.startSeqId ||
        block.startSeqId >= static_cast<int32_t>(seqDBIndexCache->seqLines.size()) ||
        block.endSeqId > static_cast<int32_t>(seqDBIndexCache->seqLines.size())) {
        std::ostringstream oss;
        oss << "The SeedDB index cache is corrupt. The block's startSeqId or endSeqId "
            << "are not valid in SeedDBIndexCache. "
            << "blockId = " << blockId << ", startSeqId = " << block.startSeqId
            << ", endSeqId = " << block.endSeqId;
        throw std::runtime_error(oss.str());
    }

    // Create a vector of sequence IDs which will be collected.
    std::vector<int32_t> seqIdsToFetch(block.endSeqId - block.startSeqId);
    std::iota(seqIdsToFetch.begin(), seqIdsToFetch.end(), block.startSeqId);

    return GetSeqDBContiguousParts(seqDBIndexCache, seqIdsToFetch);
}

std::vector<ContiguousFilePart> GetSeqDBContiguousParts(
    const std::shared_ptr<PacBio::Pancake::SeqDBIndexCache>& seqDBIndexCache,
    const std::vector<std::string>& seqNamesToFetch)
{
    HeaderLookupType headerToOrdinalId;
    ComputeSeqDBIndexHeaderLookup(*seqDBIndexCache, headerToOrdinalId);

    std::vector<std::int32_t> seqIdsToFetch;

    for (const auto& seqName : seqNamesToFetch) {
        auto it = headerToOrdinalId.find(seqName);
        if (it == headerToOrdinalId.end()) {
            throw std::runtime_error("(GetSeqDBContiguousParts) Cannot find seq name '" + seqName +
                                     "' in the provided seqDBIndexCache.");
        }
        auto id = it->second;
        seqIdsToFetch.emplace_back(id);
    }

    return GetSeqDBContiguousParts(seqDBIndexCache, seqIdsToFetch);
}

std::vector<ContiguousFilePart> GetSeqDBContiguousParts(
    const std::shared_ptr<PacBio::Pancake::SeqDBIndexCache>& seqDBIndexCache,
    std::vector<int32_t> seqIdsToFetch)
{
    // Sort the sequences by their offset in the input File.
    std::sort(seqIdsToFetch.begin(), seqIdsToFetch.end(),
              [&seqDBIndexCache](const auto& a, const auto& b) {
                  return seqDBIndexCache->GetSeqLine(a).fileOffset <
                         seqDBIndexCache->GetSeqLine(b).fileOffset;
              });

    // Sequences in the block might not be stored contiguously in the file,
    // for example if a user has permuted or filtered the DB.
    // We will collect all contiguous stretches of bytes here, and then
    // fetch those parts later.
    std::vector<ContiguousFilePart> contiguousParts;

    auto AddContiguousPart = [&](const SeqDBSequenceLine& sl) {
        contiguousParts.emplace_back(
            ContiguousFilePart{sl.fileId, sl.fileOffset, sl.fileOffset + sl.numBytes, {sl.seqId}});
    };

    for (const auto& ordId : seqIdsToFetch) {
        const auto& sl = seqDBIndexCache->GetSeqLine(ordId);

        if (contiguousParts.empty()) {
            AddContiguousPart(sl);

        } else if (sl.fileId != contiguousParts.back().fileId) {
            AddContiguousPart(sl);

        } else if (sl.fileOffset == contiguousParts.back().endOffset) {
            contiguousParts.back().endOffset += sl.numBytes;
            contiguousParts.back().seqIds.emplace_back(sl.seqId);

        } else if (sl.fileOffset > contiguousParts.back().endOffset ||
                   (sl.fileOffset + sl.numBytes) <= contiguousParts.back().startOffset) {
            // Allow out of order byte spans, as long as there is no overlap.
            AddContiguousPart(sl);

        } else {
            // An overlap occurred.
            const auto& last = contiguousParts.back();
            std::ostringstream oss;
            oss << "Invalid SeqLine object in the block, overlapping other SeqLine objects in "
                   "terms of the file offset. Last ContiguousFilePart span: {"
                << last.fileId << ", " << last.startOffset << ", " << last.endOffset
                << "}, last SeqLine: {" << sl.seqId << ", " << sl.header << ", " << sl.fileId
                << ", " << sl.fileOffset << ", " << sl.numBytes << ", " << sl.numBases << "}.";
            throw std::runtime_error(oss.str());
        }
    }

    return contiguousParts;
}

std::ostream& operator<<(std::ostream& os, const PacBio::Pancake::SeqDBIndexCache& r)
{
    os << "V\t" << r.version << "\n";
    os << "C\t" << r.compressionLevel << "\n";
    for (const auto& fl : r.fileLines) {
        os << "F"
           << "\t" << fl.fileId << "\t" << fl.filename << "\t" << fl.numSequences << "\t"
           << fl.numBytes << "\t" << fl.numCompressedBases << "\n";
    }
    for (const auto& sl : r.seqLines) {
        os << "S"
           << "\t" << sl.seqId << "\t" << sl.header << "\t" << sl.fileId << "\t" << sl.fileOffset
           << "\t" << sl.numBytes << "\t" << sl.numBases << "\t" << sl.ranges.size();
        for (const auto& r : sl.ranges) {
            os << "\t" << r.start << "\t" << r.end;
        }
        os << "\n";
    }
    for (const auto& bl : r.blockLines) {
        os << "B"
           << "\t" << bl.blockId << "\t" << bl.startSeqId << "\t" << bl.endSeqId << "\t"
           << bl.numBytes << "\t" << bl.numBases << "\n";
    }
    return os;
}

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

std::unique_ptr<PacBio::Pancake::SeqDBIndexCache> FilterSeqDBIndexCache(
    const SeqDBIndexCache& inSeqDBCache, const SamplingType& samplingType,
    const int64_t sampledBases, const int64_t randomSeed, const FilterListType& filterType,
    const std::unordered_set<std::string>& filterList, const bool doNormalization,
    const int32_t normBlockSize, const std::string& outIndexFilename)
{
    std::unique_ptr<PacBio::Pancake::SeqDBIndexCache> filteredSeqDBCache =
        std::make_unique<PacBio::Pancake::SeqDBIndexCache>();

    // Set the new filename to be the same as the old one.
    if (outIndexFilename.empty()) {
        filteredSeqDBCache->indexFilename = inSeqDBCache.indexFilename;
    } else {
        filteredSeqDBCache->indexFilename = outIndexFilename;
    }

    // Initialize the file information, version and compression level.
    SplitPath(filteredSeqDBCache->indexFilename, filteredSeqDBCache->indexParentFolder,
              filteredSeqDBCache->indexBasename);
    filteredSeqDBCache->version = inSeqDBCache.version;
    filteredSeqDBCache->compressionLevel = inSeqDBCache.compressionLevel;

    // The data will not be copied (only the index), so the file lines are the same.
    filteredSeqDBCache->fileLines = inSeqDBCache.fileLines;

    // Filter the sequence lines.
    PerformSeqDBSequenceLineSampling(filteredSeqDBCache->seqLines, inSeqDBCache.seqLines,
                                     samplingType, sampledBases, randomSeed, filterList,
                                     filterType);

    if (doNormalization) {
        NormalizeSeqDBIndexCache(*filteredSeqDBCache, normBlockSize);
    }

    return filteredSeqDBCache;
}

}  // namespace Pancake
}  // namespace PacBio
