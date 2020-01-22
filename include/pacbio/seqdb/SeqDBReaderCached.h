// Author: Ivan Sovic

#ifndef PANCAKE_SEQDB_READER_CACHED_H
#define PANCAKE_SEQDB_READER_CACHED_H

#include <pacbio/seqdb/FastaSequenceId.h>
#include <pacbio/seqdb/SeqDBIndexCache.h>
#include <pacbio/seqdb/SeqDBReaderCached.h>
#include <memory>
#include <string>

namespace PacBio {
namespace Pancake {

class SeqDBReaderCached
{
public:
    SeqDBReaderCached(std::shared_ptr<PacBio::Pancake::SeqDBIndexCache>& seedDBCache,
                      int32_t blockId);
    ~SeqDBReaderCached();

    const FastaSequenceId& GetSequence(int32_t seqId) const;
    const FastaSequenceId& GetSequence(const std::string& seqName) const;
    const std::vector<PacBio::Pancake::FastaSequenceId>& records() const { return records_; }

private:
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBIndexCache_;
    int32_t blockId_;
    std::vector<PacBio::Pancake::FastaSequenceId> records_;

    // Info to allow random access.
    std::unordered_map<std::string, int32_t> headerToOrdinalId_;
    std::unordered_map<int32_t, int32_t> seqIdToOrdinalId_;
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_SEQDB_READER_CACHED_H