// Authors: Ivan Sovic

#include <pacbio/seqdb/SeqDBIndexCache.h>
#include <pacbio/seqdb/Util.h>
#include <limits>
#include <sstream>

namespace PacBio {
namespace Pancake {

std::unique_ptr<PacBio::Pancake::SeqDBIndexCache> LoadSeqDBIndexCache(
    const std::string& indexFilename)
{
    std::ifstream ifs(indexFilename);
    return LoadSeqDBIndexCache(ifs, indexFilename);
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
