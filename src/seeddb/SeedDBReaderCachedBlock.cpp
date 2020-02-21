// Authors: Ivan Sovic

#include <pacbio/seeddb/SeedDBReader.h>
#include <pacbio/seeddb/SeedDBReaderCachedBlock.h>
#include <iostream>
#include <sstream>

namespace PacBio {
namespace Pancake {

SeedDBReaderCachedBlock::SeedDBReaderCachedBlock(
    std::shared_ptr<PacBio::Pancake::SeedDBIndexCache>& seedDBCache)
    : indexCache_(seedDBCache), blockId_(-1)
{
    // Sanity check.
    if (indexCache_->fileLines.empty())
        throw std::runtime_error("There are no file specifications in the input index file.");
    if (indexCache_->seedLines.empty())
        throw std::runtime_error("There are no sequences in the input index file.");
    if (indexCache_->blockLines.empty())
        throw std::runtime_error("There are no blocks in the input index file.");
}

SeedDBReaderCachedBlock::SeedDBReaderCachedBlock(
    std::shared_ptr<PacBio::Pancake::SeedDBIndexCache>& seedDBCache, int32_t blockId)
    : indexCache_(seedDBCache), blockId_(blockId)
{
    // Sanity check.
    if (indexCache_->fileLines.empty())
        throw std::runtime_error("There are no file specifications in the input index file.");
    if (indexCache_->seedLines.empty())
        throw std::runtime_error("There are no sequences in the input index file.");
    if (indexCache_->blockLines.empty())
        throw std::runtime_error("There are no blocks in the input index file.");

    LoadBlock(blockId);
}

SeedDBReaderCachedBlock::~SeedDBReaderCachedBlock() = default;

void SeedDBReaderCachedBlock::LoadBlock(int32_t blockId)
{
    blockId_ = blockId;

    // Form the contiguous blocks for loading.
    std::vector<ContiguousFilePart> parts = GetSeedDBContiguousParts(indexCache_, blockId_);
    const auto& bl = indexCache_->GetBlockLine(blockId_);

    // Preallocate the data size.
    int64_t totalSize = 0;
    for (int32_t sId = bl.startSeqId; sId < bl.endSeqId; ++sId) {
        const auto& sl = indexCache_->GetSeedsLine(sId);
        totalSize += sl.numSeeds;
    }
    data_.resize(totalSize);
    records_.resize(bl.endSeqId - bl.startSeqId);

    int64_t currDataPos = 0;
    int64_t currRecord = 0;
    for (const auto& part : parts) {
        // Open the file and position to the correct offset.
        const auto& fl = indexCache_->GetFileLine(part.fileId);
        const std::string actualPath = JoinPath(indexCache_->indexParentFolder, fl.filename);
        std::unique_ptr<FILE, FileDeleter> fp = PacBio::Pancake::OpenFile(actualPath.c_str(), "rb");
        const int32_t rv = fseek(fp.get(), part.startOffset, SEEK_SET);
        if (rv) {
            throw std::runtime_error("Could not fseek to position: " +
                                     std::to_string(part.startOffset));
        }

        // Load the bytes.
        const int64_t itemsToRead = (part.endOffset - part.startOffset) / 16;
        const int64_t numItemsRead =
            fread(&data_[currDataPos], sizeof(__int128), itemsToRead, fp.get());
        if (itemsToRead != numItemsRead) {
            std::ostringstream oss;
            oss << "(SeqDBReaderCachedBlock) Could not read data for the following part: "
                << "fileId = " << part.fileId << ", offsetStart = " << part.startOffset
                << ", offsetEnd = " << part.endOffset << ", startId = " << part.startId
                << ", endId = " << part.endId;
            throw std::runtime_error(oss.str());
        }

        // Create the records, and add them to the lookups.
        int64_t seqStart = currDataPos;
        for (int32_t id = part.startId; id < part.endId; ++id, ++currRecord) {
            const auto& sl = indexCache_->GetSeedsLine(id);
            records_[currRecord] =
                SequenceSeedsCached(sl.header, &data_[seqStart], sl.numSeeds, sl.seqId);
            headerToOrdinalId_[sl.header] = currRecord;
            seqIdToOrdinalId_[sl.seqId] = currRecord;
            seqStart += sl.numSeeds;
        }

        // Increment the storage location for the next part.
        currDataPos += numItemsRead;
    }
}

const SequenceSeedsCached& SeedDBReaderCachedBlock::GetSeedsForSequence(int32_t seqId) const
{
    auto it = seqIdToOrdinalId_.find(seqId);
    if (it == seqIdToOrdinalId_.end()) {
        std::ostringstream oss;
        oss << "(SeedDBReaderCachedBlock) Invalid seqId, not found in block " << blockId_
            << ". seqId = " << seqId << ", records_.size() = " << records_.size();
        throw std::runtime_error(oss.str());
    }
    int32_t ordinalId = it->second;
    return records_[ordinalId];
}

const SequenceSeedsCached& SeedDBReaderCachedBlock::GetSeedsForSequence(
    const std::string& seqName) const
{
    auto it = headerToOrdinalId_.find(seqName);
    if (it == headerToOrdinalId_.end()) {
        std::ostringstream oss;
        oss << "(SeedDBReaderCachedBlock) Invalid seqName, not found in block " << blockId_
            << ". seqName = '" << seqName << ".";
        throw std::runtime_error(oss.str());
    }
    int32_t ordinalId = it->second;
    return records_[ordinalId];
}

}  // namespace Pancake
}  // namespace PacBio
