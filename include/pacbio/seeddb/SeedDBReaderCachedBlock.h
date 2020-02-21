// Author: Ivan Sovic

#ifndef PANCAKE_SEEDDB_READER_CACHED_BLOCK_H
#define PANCAKE_SEEDDB_READER_CACHED_BLOCK_H

#include <pacbio/seeddb/SeedDBIndexCache.h>
#include <pacbio/seeddb/SeedDBReaderCachedBlock.h>
#include <pacbio/seeddb/SequenceSeedsCached.h>
#include <pacbio/seqdb/Util.h>
#include <memory>
#include <ostream>
#include <string>

namespace PacBio {
namespace Pancake {

class SeedDBReaderCachedBlock
{
public:
    SeedDBReaderCachedBlock(std::shared_ptr<PacBio::Pancake::SeedDBIndexCache>& seedDBCache);
    SeedDBReaderCachedBlock(std::shared_ptr<PacBio::Pancake::SeedDBIndexCache>& seedDBCache,
                            int32_t blockId);
    ~SeedDBReaderCachedBlock();

    void LoadBlock(int32_t blockId);
    const SequenceSeedsCached& GetSeedsForSequence(int32_t seqId) const;
    const SequenceSeedsCached& GetSeedsForSequence(const std::string& seqName) const;
    const std::vector<PacBio::Pancake::SequenceSeedsCached>& records() const { return records_; }

private:
    std::shared_ptr<PacBio::Pancake::SeedDBIndexCache> indexCache_;
    int32_t blockId_;
    std::vector<__int128> data_;
    std::vector<PacBio::Pancake::SequenceSeedsCached> records_;

    // Info to allow random access.
    std::unordered_map<std::string, int32_t> headerToOrdinalId_;
    std::unordered_map<int32_t, int32_t> seqIdToOrdinalId_;
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_SEEDDB_READER_CACHED_BLOCK_H