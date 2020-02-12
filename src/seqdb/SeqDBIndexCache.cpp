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

    while ((numRead = getline(&line, &lineLen, fpIn)) != -1) {
        if (lineLen <= 0) {
            Cleanup();
            continue;
        }
        SeqDBFileLine fl;
        SeqDBSequenceLine sl;
        SeqDBBlockLine bl;
        int32_t numRanges = 0;
        int32_t ordinalId = 0;
        int32_t n = 0;
        int32_t readOffset = 0;
        int32_t offset = 0;
        const char token = line[0];
        switch (token) {
            case 'V':
                n = sscanf(&line[1], "%s", buff);
                cache->version = buff;
                if (n != 1)
                    throw std::runtime_error("Problem parsing line: '" + std::string(line) + "'.");
                break;
            case 'C':
                n = sscanf(&line[1], "%d", &(cache->compressionLevel));
                if (n != 1)
                    throw std::runtime_error("Problem parsing line: '" + std::string(line) + "'.");
                break;
            case 'F':
                n = sscanf(&line[1], "%d %s %d %lld %lld", &(fl.fileId), buff, &(fl.numSequences),
                           &(fl.numBytes), &(fl.numCompressedBases));
                if (n != 5)
                    throw std::runtime_error("Problem parsing line: '" + std::string(line) + "'.");
                fl.filename = buff;
                cache->fileLines.emplace_back(fl);
                break;
            case 'S':
                n = sscanf(&line[1], "%d %s %d %lld %d %d %d%n", &(sl.seqId), buff, &(sl.fileId),
                           &(sl.fileOffset), &(sl.numBytes), &(sl.numBases), &(numRanges),
                           &readOffset);
                if (n != 7)
                    throw std::runtime_error("Problem parsing line: '" + std::string(line) + "'.");
                sl.header = buff;
                offset = readOffset + 1;
                for (int32_t i = 0; i < numRanges; ++i) {
                    Range r;
                    n = sscanf(&line[offset], "%d %d%n", &r.start, &r.end, &readOffset);
                    if (n != 2)
                        throw std::runtime_error("Problem parsing line: '" + std::string(line) +
                                                 "'.");
                    offset += readOffset + 1;
                    sl.ranges.emplace_back(r);
                }
                ordinalId = cache->seqLines.size();
                cache->seqLines.emplace_back(sl);
                cache->headerToOrdinalId[sl.header] = ordinalId;
                cache->seqIdToOrdinalId[sl.seqId] = ordinalId;
                break;
            case 'B':
                n = sscanf(&line[1], "%d %d %d %lld %lld", &(bl.blockId), &(bl.startSeqId),
                           &(bl.endSeqId), &(bl.numBytes), &(bl.numBases));
                if (n != 5)
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
        int32_t ordinalId = 0;

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
                ordinalId = cache->seqLines.size();
                cache->seqLines.emplace_back(sl);
                // Add the new sequence to the lookups.
                cache->headerToOrdinalId[sl.header] = ordinalId;
                cache->seqIdToOrdinalId[sl.seqId] = ordinalId;
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
