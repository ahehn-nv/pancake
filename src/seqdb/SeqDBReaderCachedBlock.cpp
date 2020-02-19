// Authors: Ivan Sovic

#include <pacbio/seqdb/SeqDBReader.h>
#include <pacbio/seqdb/SeqDBReaderCachedBlock.h>
#include <pacbio/seqdb/Twobit.h>
#include <pacbio/seqdb/Util.h>
#include <sstream>

namespace PacBio {
namespace Pancake {

SeqDBReaderCachedBlock::SeqDBReaderCachedBlock(
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache>& seqDBCache)
    : seqDBIndexCache_(seqDBCache)
{
    // Sanity check.
    if (seqDBIndexCache_->fileLines.empty())
        throw std::runtime_error("There are no file specifications in the input index file.");
    if (seqDBIndexCache_->seqLines.empty())
        throw std::runtime_error("There are no sequences in the input index file.");
    if (seqDBIndexCache_->blockLines.empty())
        throw std::runtime_error("There are no blocks in the input index file.");
}

SeqDBReaderCachedBlock::SeqDBReaderCachedBlock(
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache>& seqDBCache, int32_t blockId)
    : seqDBIndexCache_(seqDBCache), blockId_(blockId)
{
    // Sanity check.
    if (seqDBIndexCache_->fileLines.empty())
        throw std::runtime_error("There are no file specifications in the input index file.");
    if (seqDBIndexCache_->seqLines.empty())
        throw std::runtime_error("There are no sequences in the input index file.");
    if (seqDBIndexCache_->blockLines.empty())
        throw std::runtime_error("There are no blocks in the input index file.");

    LoadBlock(blockId_);
}

SeqDBReaderCachedBlock::~SeqDBReaderCachedBlock() = default;

void SeqDBReaderCachedBlock::LoadBlock(int32_t blockId)
{
    if (seqDBIndexCache_->compressionLevel == 0) {
        return LoadBlockUncompressed_(blockId);
    }
    return LoadBlockCompressed_(blockId);
}

void SeqDBReaderCachedBlock::LoadBlockCompressed_(int32_t blockId)
{
    blockId_ = blockId;

    // Form the contiguous blocks for loading.
    std::vector<ContiguousFilePart> parts = GetSeqDBContiguousParts(seqDBIndexCache_, blockId_);
    const auto& bl = seqDBIndexCache_->GetBlockLine(blockId_);

    // Preallocate the data size.
    int64_t totalSize = 0;
    for (int32_t sId = bl.startSeqId; sId < bl.endSeqId; ++sId) {
        const auto& sl = seqDBIndexCache_->GetSeqLine(sId);
        totalSize += sl.numBases;
    }
    data_.resize(totalSize);
    records_.resize(bl.endSeqId - bl.startSeqId);

    int64_t currDataPos = 0;
    int64_t currRecord = 0;
    for (const auto& part : parts) {
        // Open the file and position to the correct offset.
        const auto& fl = seqDBIndexCache_->GetFileLine(part.fileId);
        const std::string actualPath = JoinPath(seqDBIndexCache_->indexParentFolder, fl.filename);
        std::unique_ptr<FILE, FileDeleter> fp = PacBio::Pancake::OpenFile(actualPath.c_str(), "rb");
        const int32_t rv = fseek(fp.get(), part.startOffset, SEEK_SET);
        if (rv) {
            throw std::runtime_error("Could not fseek to position: " +
                                     std::to_string(part.startOffset));
        }

        // Load the bytes.
        const int64_t itemsToRead = (part.endOffset - part.startOffset);
        std::vector<uint8_t> tempData;
        tempData.resize(itemsToRead);
        const int64_t numItemsRead = fread(&tempData[0], sizeof(uint8_t), itemsToRead, fp.get());
        if (itemsToRead != numItemsRead) {
            std::ostringstream oss;
            oss << "(SeqDBReaderCachedBlock) Could not read data for the following part: "
                << "fileId = " << part.fileId << ", offsetStart = " << part.startOffset
                << ", offsetEnd = " << part.endOffset << ", startId = " << part.startId
                << ", endId = " << part.endId;
            throw std::runtime_error(oss.str());
        }

        // Position of a current record in the data_ vector.
        int64_t seqStart = currDataPos;

        // File offset of the first record in the part. Needed to normalize
        // the access to the tempData vector.
        int64_t startByteOffset = part.startOffset;

        // Decompress and create the records, and add them to the lookups.
        for (int32_t id = part.startId; id < part.endId; ++id, ++currRecord) {
            const auto& sl = seqDBIndexCache_->GetSeqLine(id);
            const int64_t firstByte = sl.fileOffset - startByteOffset;
            DecompressSequence(&tempData[firstByte], tempData.size(), sl.numBases, sl.ranges,
                               &data_[seqStart]);
            records_[currRecord] = FastaSequenceCached{
                sl.header, reinterpret_cast<const char*>(&data_[seqStart]), sl.numBases, sl.seqId};
            headerToOrdinalId_[sl.header] = currRecord;
            seqIdToOrdinalId_[sl.seqId] = currRecord;
            seqStart += sl.numBases;
        }

        // Increment the storage location for the next part.
        currDataPos += numItemsRead;
    }
}

void SeqDBReaderCachedBlock::LoadBlockUncompressed_(int32_t blockId)
{
    blockId_ = blockId;

    // Form the contiguous blocks for loading.
    std::vector<ContiguousFilePart> parts = GetSeqDBContiguousParts(seqDBIndexCache_, blockId_);
    const auto& bl = seqDBIndexCache_->GetBlockLine(blockId_);

    // Preallocate the data size.
    // int64_t totalSize = 0;
    // for (const auto& part : parts) {
    //     totalSize += (part.endOffset - part.startOffset);
    // }
    int64_t totalSize = 0;
    for (int32_t sId = bl.startSeqId; sId < bl.endSeqId; ++sId) {
        const auto& sl = seqDBIndexCache_->GetSeqLine(sId);
        totalSize += sl.numBases;
    }
    data_.resize(totalSize);
    records_.resize(bl.endSeqId - bl.startSeqId);

    int64_t currDataPos = 0;
    int64_t currRecord = 0;
    for (const auto& part : parts) {
        // Open the file and position to the correct offset.
        const auto& fl = seqDBIndexCache_->GetFileLine(part.fileId);
        const std::string actualPath = JoinPath(seqDBIndexCache_->indexParentFolder, fl.filename);
        std::unique_ptr<FILE, FileDeleter> fp = PacBio::Pancake::OpenFile(actualPath.c_str(), "rb");
        const int32_t rv = fseek(fp.get(), part.startOffset, SEEK_SET);
        if (rv) {
            throw std::runtime_error("Could not fseek to position: " +
                                     std::to_string(part.startOffset));
        }

        // Load the bytes.
        const int64_t itemsToRead = (part.endOffset - part.startOffset);
        const int64_t numItemsRead =
            fread(&data_[currDataPos], sizeof(uint8_t), itemsToRead, fp.get());
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
            const auto& sl = seqDBIndexCache_->GetSeqLine(id);
            records_[currRecord] = FastaSequenceCached{
                sl.header, reinterpret_cast<const char*>(&data_[seqStart]), sl.numBases, sl.seqId};
            headerToOrdinalId_[sl.header] = currRecord;
            seqIdToOrdinalId_[sl.seqId] = currRecord;
            seqStart += sl.numBases;
        }

        // Increment the storage location for the next part.
        currDataPos += numItemsRead;
    }
}

const FastaSequenceCached& SeqDBReaderCachedBlock::GetSequence(int32_t seqId) const
{
    auto it = seqIdToOrdinalId_.find(seqId);
    if (it == seqIdToOrdinalId_.end()) {
        std::ostringstream oss;
        oss << "(SeqDBReaderCachedBlock) Invalid seqId, not found in block " << blockId_
            << ". seqId = " << seqId << ", records_.size() = " << records_.size();
        throw std::runtime_error(oss.str());
    }
    int32_t ordinalId = it->second;
    return records_[ordinalId];
}

const FastaSequenceCached& SeqDBReaderCachedBlock::GetSequence(const std::string& seqName) const
{
    auto it = headerToOrdinalId_.find(seqName);
    if (it == headerToOrdinalId_.end()) {
        std::ostringstream oss;
        oss << "(SeqDBReaderCachedBlock) Invalid seqName, not found in block " << blockId_
            << ". seqName = '" << seqName << ".";
        throw std::runtime_error(oss.str());
    }
    int32_t ordinalId = it->second;
    return records_[ordinalId];
}

}  // namespace Pancake
}  // namespace PacBio
