// Author: Ivan Sovic

#ifndef PANCAKE_SEQDB_READER_H
#define PANCAKE_SEQDB_READER_H

#include <pacbio/seqdb/SeqDBReaderBase.h>
#include <pacbio/seqdb/SeqDBWriterBase.h>
#include <pacbio/seqdb/Util.h>
#include <pbbam/FastaSequence.h>
#include <memory>
#include <ostream>
#include <string>

namespace PacBio {
namespace Pancake {

class SeqDBReader : SeqDBReaderBase
{
public:
    SeqDBReader(std::shared_ptr<PacBio::Pancake::SeqDBIndexCache>& seqDBCache);
    ~SeqDBReader();
    bool GetSequence(BAM::FastaSequence& record, int64_t seqId) override;
    bool GetSequence(BAM::FastaSequence& record, const std::string& seqName) override;
    bool GetNext(BAM::FastaSequence& record) override;
    bool GetNextBatch(std::vector<BAM::FastaSequence>& records, int64_t batchSize) override;
    bool GetBlock(std::vector<BAM::FastaSequence>& records, int32_t blockId) override;

    bool JumpTo(int64_t seqId) override;
    bool JumpTo(const std::string& seqName) override;

private:
    using FilePtr = std::unique_ptr<FILE, FileDeleter>;
    class OpenFileHandler
    {
    public:
        FilePtr fp = nullptr;
        int32_t fileId = -1;
        int64_t pos = 0;
        int32_t nextOrdinalId = 0;
    };
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBIndexCache_;
    OpenFileHandler fileHandler_;

    void AccessLocation_(OpenFileHandler& fileHandler,
                         const std::vector<PacBio::Pancake::SeqDBFileLine>& fileLines,
                         const std::string& indexParentFolder, int32_t fileId,
                         int32_t nextOrdinalId, int64_t offset) const;

    void LoadAndUnpackSequence_(BAM::FastaSequence& record, OpenFileHandler& fileHandler,
                                const std::vector<PacBio::Pancake::SeqDBFileLine>& fileLines,
                                const std::string& indexParentFolder, const SeqDBSequenceLine& sl,
                                int32_t ordinalId, bool isCompressed) const;
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_SEQDB_READER_H