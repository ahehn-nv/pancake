// Authors: Ivan Sovic

#include <pacbio/seqdb/SeqDBReader.h>
#include <pacbio/seqdb/SeqDBReaderCached.h>
#include <sstream>

namespace PacBio {
namespace Pancake {

SeqDBReaderCached::SeqDBReaderCached(std::shared_ptr<PacBio::Pancake::SeqDBIndexCache>& seqDBCache,
                                     int32_t blockId)
    : seqDBIndexCache_(seqDBCache), blockId_(blockId)
{
    // Sanity check.
    if (seqDBIndexCache_->fileLines.empty())
        throw std::runtime_error("There are no file specifications in the input index file.");
    if (seqDBIndexCache_->seqLines.empty())
        throw std::runtime_error("There are no sequences in the input index file.");
    if (seqDBIndexCache_->blockLines.empty())
        throw std::runtime_error("There are no blocks in the input index file.");

    // Create a SeedDB reader.
    PacBio::Pancake::SeqDBReader reader(seqDBCache);

    // Fetch the block of data.
    reader.GetBlock(records_, blockId);

    for (int32_t i = 0; i < static_cast<int32_t>(records_.size()); ++i) {
        headerToOrdinalId_[records_[i].Name()] = i;
        seqIdToOrdinalId_[records_[i].Id()] = i;
    }
}

SeqDBReaderCached::~SeqDBReaderCached() = default;

const FastaSequenceId& SeqDBReaderCached::GetSequence(int32_t seqId) const
{
    auto it = seqIdToOrdinalId_.find(seqId);
    if (it == seqIdToOrdinalId_.end()) {
        std::ostringstream oss;
        oss << "(SeqDBReaderCached) Invalid seqId, not found in block " << blockId_
            << ". seqId = " << seqId << ", records_.size() = " << records_.size();
        throw std::runtime_error(oss.str());
    }
    int32_t ordinalId = it->second;
    return records_[ordinalId];
}

const FastaSequenceId& SeqDBReaderCached::GetSequence(const std::string& seqName) const
{
    auto it = headerToOrdinalId_.find(seqName);
    if (it == headerToOrdinalId_.end()) {
        std::ostringstream oss;
        oss << "(SeqDBReaderCached) Invalid seqName, not found in block " << blockId_
            << ". seqName = '" << seqName << ".";
        throw std::runtime_error(oss.str());
    }
    int32_t ordinalId = it->second;
    return records_[ordinalId];
}

}  // namespace Pancake
}  // namespace PacBio
