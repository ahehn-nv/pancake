// Author: Ivan Sovic

#ifndef PANCAKE_SEEDDB_INDEX_CACHE_H
#define PANCAKE_SEEDDB_INDEX_CACHE_H

#include <pacbio/seeddb/SeedDBParameters.h>
#include <pacbio/seqdb/Range.h>
#include <cstdint>
#include <fstream>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace PacBio {
namespace Pancake {

/// \brief      Container, describes a seeds file which accompanies the SeedDB index.
///
class SeedDBFileLine
{
public:
    int32_t fileId = 0;
    std::string filename;
    int32_t numSequences = 0;
    int64_t numBytes = 0;
};

/// \brief      Container, index information for a particular sequence's set of seeds.
///
class SeedDBSeedsLine
{
public:
    int32_t seqId = 0;
    std::string header;
    int32_t fileId = 0;
    int64_t fileOffset = 0;
    int64_t numBytes = 0;
    int32_t numBases = 0;
    int32_t numSeeds = 0;
};

/// \brief      Container, index information for a particular block of sequences.
///
class SeedDBBlockLine
{
public:
    int32_t blockId = 0;
    int32_t startSeqId = -1;
    int32_t endSeqId = -1;
    int64_t numBytes = 0;

    int32_t Span() const { return endSeqId - startSeqId; }
};

/// \brief      Container used to store all the indexing information from the DB.
///
class SeedDBIndexCache
{
public:
    std::string indexFilename;
    std::string indexParentFolder;
    std::string indexBasename;
    std::string version{"unknown"};
    PacBio::Pancake::SeedDB::SeedDBParameters seedParams;
    std::vector<SeedDBFileLine> fileLines;
    std::vector<SeedDBSeedsLine> seedLines;
    std::vector<SeedDBBlockLine> blockLines;
    // Header to ordinal ID in the seedLines;
    std::unordered_map<std::string, int32_t> headerToOrdinalId;
    // In case the cache represents a sliced portion of the index, the
    // sequence ID can be different than the order of appearance in the seedLines
    // vector. This lookup relates the seqID to the ordinal ID.
    std::unordered_map<int32_t, int32_t> seqIdToOrdinalId;

    SeedDBIndexCache() = default;
    ~SeedDBIndexCache() = default;
    friend std::ostream& operator<<(std::ostream& os, const PacBio::Pancake::SeedDBIndexCache& r);
};

/// \brief Loads the SeedDB index from file.
///
/// \param[in]  indexFilename           Path to the .seqdb index file.
///
/// \returns Pointer to the parsed index.
///
std::unique_ptr<PacBio::Pancake::SeedDBIndexCache> LoadSeedDBIndexCache(
    const std::string& indexFilename);

std::unique_ptr<PacBio::Pancake::SeedDBIndexCache> LoadSeedDBIndexCache(
    std::istream& is, const std::string& indexFilename);

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_SEQDB_CACHE_H