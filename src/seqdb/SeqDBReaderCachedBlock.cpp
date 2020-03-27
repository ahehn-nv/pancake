// Authors: Ivan Sovic

#include <pacbio/seqdb/SeqDBReader.h>
#include <pacbio/seqdb/SeqDBReaderCachedBlock.h>
#include <pacbio/seqdb/Twobit.h>
#include <pacbio/seqdb/Util.h>
#include <iostream>
#include <sstream>

namespace PacBio {
namespace Pancake {

SeqDBReaderCachedBlock::SeqDBReaderCachedBlock(
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache>& seqDBCache)
    : seqDBIndexCache_(seqDBCache)
{
    ValidateSeqDBIndexCache(seqDBCache);
}

SeqDBReaderCachedBlock::~SeqDBReaderCachedBlock() = default;

void SeqDBReaderCachedBlock::LoadBlocks(const std::vector<int32_t>& blockIds)
{
    // Form the contiguous blocks for loading.
    // First, get all the spans for the first block. These will be stored directly in the
    // accumulator called 'parts'.
    std::vector<ContiguousFilePart> parts = GetSeqDBContiguousParts(seqDBIndexCache_, blockIds[0]);
    if (parts.empty()) {
        std::runtime_error(
            "Unknown error occurred, the GetSeqDBContiguousParts returned empty in "
            "LoadBlockUncompressed_.");
    }
    // Next, get all other parts for all other blocks. If possible, extend
    // the previously added parts, otherwise just append.
    for (size_t i = 1; i < blockIds.size(); ++i) {
        std::vector<ContiguousFilePart> newParts =
            GetSeqDBContiguousParts(seqDBIndexCache_, blockIds[i]);

        for (const auto& newPart : newParts) {
            if (newPart.CanAppendTo(parts.back())) {
                parts.back().ExtendWith(newPart);
            } else {
                parts.emplace_back(newPart);
            }
        }
    }

    // Preallocate the data size.
    int64_t totalSize = 0;
    for (const auto& blockId : blockIds) {
        const auto& bl = seqDBIndexCache_->GetBlockLine(blockId);
        for (int32_t sId = bl.startSeqId; sId < bl.endSeqId; ++sId) {
            const auto& sl = seqDBIndexCache_->GetSeqLine(sId);
            totalSize += sl.numBases;
        }
    }
    data_.resize(totalSize);

    // Preallocate the space for all the records.
    const auto& blFirst = seqDBIndexCache_->GetBlockLine(blockIds.front());
    const auto& blLast = seqDBIndexCache_->GetBlockLine(blockIds.back());
    records_.resize(blLast.endSeqId - blFirst.startSeqId);

    // Actually load the data.
    if (seqDBIndexCache_->compressionLevel == 0) {
        return LoadBlockUncompressed_(parts);
    }
    return LoadBlockCompressed_(parts);
}

void SeqDBReaderCachedBlock::LoadBlockCompressed_(const std::vector<ContiguousFilePart>& parts)
{
    // Position of a current record in the data_ vector.
    int64_t seqStart = 0;
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
                << ", endId = " << part.endId << ", itemsToRead = " << itemsToRead
                << ", numItemsRead = " << numItemsRead;
            throw std::runtime_error(oss.str());
        }

        // File offset of the first record in the part. Needed to normalize
        // the access to the tempData vector.
        int64_t startByteOffset = part.startOffset;

        // Decompress and create the records, and add them to the lookups.
        for (int32_t id = part.startId; id < part.endId; ++id, ++currRecord) {
            const auto& sl = seqDBIndexCache_->GetSeqLine(id);
            const int64_t firstByte = sl.fileOffset - startByteOffset;
            DecompressSequence(&tempData[firstByte], sl.numBytes, sl.numBases, sl.ranges,
                               &data_[seqStart]);
            records_[currRecord] = FastaSequenceCached{
                sl.header, reinterpret_cast<const char*>(&data_[seqStart]), sl.numBases, sl.seqId};
            headerToOrdinalId_[sl.header] = currRecord;
            seqIdToOrdinalId_[sl.seqId] = currRecord;
            seqStart += sl.numBases;
        }
    }
}

void SeqDBReaderCachedBlock::LoadBlockUncompressed_(const std::vector<ContiguousFilePart>& parts)
{
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
                << ", endId = " << part.endId << ", itemsToRead = " << itemsToRead
                << ", numItemsRead = " << numItemsRead;
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
        oss << "(SeqDBReaderCachedBlock) Invalid seqId, not found in the preloaded data. seqId = "
            << seqId << ", records_.size() = " << records_.size();
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
        oss << "(SeqDBReaderCachedBlock) Invalid seqName, not found in the preloaded data. seqName "
               "= '"
            << seqName << ".";
        throw std::runtime_error(oss.str());
    }
    int32_t ordinalId = it->second;
    return records_[ordinalId];
}

}  // namespace Pancake
}  // namespace PacBio
