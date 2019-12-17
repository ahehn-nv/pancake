// Author: Ivan Sovic

#ifndef PANCAKE_SEQDB_READER_BASE_H
#define PANCAKE_SEQDB_READER_BASE_H

#include <pacbio/seqdb/SeqDBIndexCache.h>
#include <pbbam/FastaSequence.h>
#include <string>
#include <vector>

namespace PacBio {
namespace Pancake {

/// \brief      Base class that can be specialized for reading a compressed or an
///             uncompressed set of sequences.
///
class SeqDBReaderBase
{
public:
    virtual ~SeqDBReaderBase() {}
    virtual bool GetSequence(BAM::FastaSequence& record, int64_t seqId) = 0;
    virtual bool GetSequence(BAM::FastaSequence& record, const std::string& seqName) = 0;
    virtual bool GetNext(BAM::FastaSequence& record) = 0;
    virtual bool GetNextBatch(std::vector<BAM::FastaSequence>& records, int64_t batchSize) = 0;
    virtual bool GetBlock(std::vector<BAM::FastaSequence>& records, int32_t blockId) = 0;
    virtual bool JumpTo(int64_t seqId) = 0;
    virtual bool JumpTo(const std::string& seqName) = 0;
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_SEQDB_READER_COMPRESSED_H