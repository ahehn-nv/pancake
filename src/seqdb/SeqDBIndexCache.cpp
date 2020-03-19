// Authors: Ivan Sovic

#include <pacbio/seqdb/SeqDBIndexCache.h>
#include <pacbio/seqdb/Util.h>
#include <array>
#include <iostream>
#include <limits>
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

    // Sequences in the block might not be stored contiguously in the file,
    // for example if a user has permuted or filtered the DB.
    // We will collect all contiguous stretches of bytes here, and then
    // fetch those parts later.
    std::vector<ContiguousFilePart> contiguousParts;

    auto AddContiguousPart = [&](const SeqDBSequenceLine& sl) {
        contiguousParts.emplace_back(ContiguousFilePart{
            sl.fileId, sl.fileOffset, sl.fileOffset + sl.numBytes, sl.seqId, sl.seqId + 1});
    };

    for (int32_t ordId = block.startSeqId; ordId < block.endSeqId; ++ordId) {
        const auto& sl = seqDBIndexCache->seqLines[ordId];

        if (contiguousParts.empty()) {
            AddContiguousPart(sl);

        } else if (sl.fileId != contiguousParts.back().fileId) {
            AddContiguousPart(sl);

        } else if (sl.fileOffset == contiguousParts.back().endOffset) {
            contiguousParts.back().endOffset += sl.numBytes;
            contiguousParts.back().endId = sl.seqId + 1;

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

}  // namespace Pancake
}  // namespace PacBio
