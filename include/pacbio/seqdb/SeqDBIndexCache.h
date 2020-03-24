// Author: Ivan Sovic

#ifndef PANCAKE_SEQDB_INDEX_CACHE_H
#define PANCAKE_SEQDB_INDEX_CACHE_H

#include <pacbio/pancake/ContiguousFilePart.h>
#include <pacbio/seqdb/Range.h>
#include <pacbio/util/CommonTypes.h>
#include <cstdint>
#include <fstream>
#include <lib/flat_hash_map/flat_hash_map.hpp>
#include <memory>
#include <random>
#include <string>
#include <unordered_map>
#include <unordered_set>
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

class SeqDBBlockLine
{
public:
    int32_t blockId = 0;
    int32_t startSeqId = -1;
    int32_t endSeqId = -1;
    int64_t numBytes = 0;
    int64_t numBases = 0;

    int32_t Span() const { return endSeqId - startSeqId; }
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
    std::vector<SeqDBBlockLine> blockLines;

    SeqDBIndexCache() = default;
    ~SeqDBIndexCache() = default;
    friend std::ostream& operator<<(std::ostream& os, const PacBio::Pancake::SeqDBIndexCache& r);

    const SeqDBSequenceLine& GetSeqLine(int32_t seqId) const;
    const SeqDBBlockLine& GetBlockLine(int32_t blockId) const;
    const SeqDBFileLine& GetFileLine(int32_t fileId) const;
};

/// \brief Writes the .seqdb file to disk.
void WriteSeqDBIndexCache(FILE* fpOut, const SeqDBIndexCache& cache);

void ComputeSeqDBIndexHeaderLookup(const PacBio::Pancake::SeqDBIndexCache& dbCache,
                                   HeaderLookupType& headerToOrdinalId);

/// \brief Given a vector of SeqDBSequenceLine objects, this function groups them in
///        creates blocks of given block size.
std::vector<SeqDBBlockLine> CreateSeqDBBlocks(const std::vector<SeqDBSequenceLine>& seqLines,
                                              int64_t blockSize);

/// \brief When a SeqDBIndexCache is subsampled or manually constructed, the seqID of each sequence
///        might not correspond to the sequence's ordinal ID. This function updates the seqID to
///        match, and also constructs the blocks of given block size.
///        It modifies the cache object in place.
void NormalizeSeqDBIndexCache(SeqDBIndexCache& cache, int64_t blockSize);

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
    const std::string& indexFilename);

std::unique_ptr<PacBio::Pancake::SeqDBIndexCache> LoadSeqDBIndexCache(
    std::istream& is, const std::string& indexFilename);

std::unique_ptr<PacBio::Pancake::SeqDBIndexCache> LoadSeqDBIndexCache(
    FILE* fpIn, const std::string& indexFilename);

std::vector<ContiguousFilePart> GetSeqDBContiguousParts(
    const std::shared_ptr<PacBio::Pancake::SeqDBIndexCache>& seqDBIndexCache, int32_t blockId);

void PerformSeqDBSequenceLineSampling(std::vector<SeqDBSequenceLine>& outSeqLines,
                                      const std::vector<SeqDBSequenceLine>& inSeqLines,
                                      const SamplingType& sampling, int64_t sampledBases,
                                      const int64_t randomSeed,
                                      const std::unordered_set<std::string>& filterList,
                                      const FilterListType& filterType);

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_SEQDB_CACHE_H