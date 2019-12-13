// Author: Ivan Sovic

#ifndef PANCAKE_SEQDB_CACHE_H
#define PANCAKE_SEQDB_CACHE_H

#include <pacbio/seqdb/Range.h>
#include <cstdint>
#include <string>

namespace PacBio {
namespace Pancake {

class SeqDBFileLine
{
public:
    int32_t fileId = 0;
    std::string filename;
    int32_t numSequences = 0;
    int64_t numBytes = 0;
    int64_t numUncompressedBases = 0;
    int64_t numCompressedBases = 0;
};

class SeqDBSequenceLine
{
public:
    int32_t seqId = 0;
    std::string header;
    int32_t numBytes = 0;
    int32_t numBases = 0;
    int32_t fileId = 0;
    int64_t fileOffset = 0;
    std::vector<PacBio::Pancake::Range> ranges;
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_SEQDB_CACHE_H