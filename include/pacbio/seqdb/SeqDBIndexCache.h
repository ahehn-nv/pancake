// Author: Ivan Sovic

#ifndef PANCAKE_SEQDB_INDEX_CACHE_H
#define PANCAKE_SEQDB_INDEX_CACHE_H

#include <pacbio/seqdb/Range.h>
#include <cstdint>
#include <fstream>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace PacBio {
namespace Pancake {

/// \brief      Container, describes a sequence file which accompanies the DB index.
///
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

/// \brief      Container, index information for a particular sequence.
///
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

/// \brief      Container used to store all the indexing information from the DB.
///
class SeqDBIndexCache
{
public:
    std::string indexFilename;
    std::string indexParentFolder;
    std::string indexBasename;
    std::string version{"unknown"};
    int32_t compressionLevel = 0;
    std::vector<SeqDBFileLine> fileLines;
    std::vector<SeqDBSequenceLine> seqLines;
    // Header to ordinal ID in the seqLines;
    std::unordered_map<std::string, int32_t> headerToOrdinalId;
    // In case the cache represents a sliced portion of the index, the
    // sequence ID can be different than the order of appearance in the seqLines
    // vector. This lookup relates the seqID to the ordinal ID.
    std::unordered_map<int32_t, int32_t> seqIdToOrdinalId;
    std::vector<PacBio::Pancake::Range> blocks;

    SeqDBIndexCache() = default;
    ~SeqDBIndexCache() = default;
    friend std::ostream& operator<<(std::ostream& os, const PacBio::Pancake::SeqDBIndexCache& r);
};

/// \brief Loads the SeqDB index from file.
///
/// \note  This function also computes the ranges of blocks give a block size
///        in bytes. If the block size is < 0, then the entire index is in
///        a single block. If block size is 0, each sequence will be in one block.

/// \param[in]  indexFilename           Path to the .seqdb index file.
/// \param[out] blockSizeInBytes        Block size.
///
/// \returns Pointer to the parsed index.
///
std::unique_ptr<PacBio::Pancake::SeqDBIndexCache> LoadSeqDBIndexCache(
    const std::string& indexFilename, int64_t blockSizeInBytes = -1);

std::unique_ptr<PacBio::Pancake::SeqDBIndexCache> LoadSeqDBIndexCache(
    std::istream& is, const std::string& indexFilename, int64_t blockSizeInBytes = -1);

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_SEQDB_CACHE_H