// Authors: Ivan Sovic

#include <pacbio/seqdb/SeqDBIndexCache.h>
#include <pacbio/seqdb/Util.h>
#include <limits>
#include <sstream>

namespace PacBio {
namespace Pancake {

std::unique_ptr<PacBio::Pancake::SeqDBIndexCache> LoadSeqDBIndexCache(
    const std::string& indexFilename, int64_t blockSizeInBytes)
{
    std::ifstream ifs(indexFilename);
    return LoadSeqDBIndexCache(ifs, indexFilename, blockSizeInBytes);
}

std::unique_ptr<PacBio::Pancake::SeqDBIndexCache> LoadSeqDBIndexCache(
    std::istream& is, const std::string& indexFilename, int64_t blockSizeInBytes)
{
    auto cache = std::make_unique<PacBio::Pancake::SeqDBIndexCache>();

    cache->indexFilename = indexFilename;
    SplitPath(indexFilename, cache->indexParentFolder, cache->indexBasename);

    // If the blockSizeInBytes is < 0, then all sequences should be in the same block.
    blockSizeInBytes =
        (blockSizeInBytes < 0) ? std::numeric_limits<int64_t>::max() : blockSizeInBytes;

    std::string line;
    char token;
    PacBio::Pancake::Range range;
    int64_t currBlockBytes = 0;
    while (std::getline(is, line)) {
        if (line.empty()) continue;

        std::istringstream iss(line);
        iss >> token;

        SeqDBFileLine fl;
        SeqDBSequenceLine sl;
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
                // Extend the block reach.
                range.end = sl.seqId + 1;
                currBlockBytes += sl.numBytes;
                // Start a new block.
                if (currBlockBytes >= blockSizeInBytes) {
                    cache->blocks.emplace_back(range);
                    range.start = sl.seqId + 1;
                    currBlockBytes = 0;
                }
                break;
            default:
                std::ostringstream oss;
                oss << "Unknown token found when parsing the index: " << token;
                throw std::runtime_error(oss.str());
                break;
        }
    }

    if (cache->seqLines.empty())
        throw std::runtime_error("There were no sequences in the input index file.");

    range.end = static_cast<int32_t>(cache->seqLines.size());
    if (range.Span() > 0) cache->blocks.emplace_back(range);

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
    return os;
}

}  // namespace Pancake
}  // namespace PacBio
