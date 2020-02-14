// Authors: Ivan Sovic

#include <pacbio/seeddb/SeedDBReader.h>
#include <pacbio/seeddb/SeedDBReaderRawBlock.h>
#include <functional>
#include <iostream>
#include <sstream>

namespace PacBio {
namespace Pancake {

SeedDBReaderRawBlock::SeedDBReaderRawBlock(
    const std::shared_ptr<PacBio::Pancake::SeedDBIndexCache>& seedDBCache)
    : seedDBIndexCache_(seedDBCache)
{
    // Sanity check.
    if (seedDBIndexCache_->fileLines.empty())
        throw std::runtime_error("There are no file specifications in the input index file.");
    if (seedDBIndexCache_->seedLines.empty())
        throw std::runtime_error("There are no sequences in the input index file.");
    if (seedDBIndexCache_->blockLines.empty())
        throw std::runtime_error("There are no blocks in the input index file.");
}

std::vector<SeedDB::SeedRaw> SeedDBReaderRawBlock::GetBlock(int32_t blockId) const
{
    // Sanity check for the sequence ID.
    if (blockId < 0 || blockId >= static_cast<int32_t>(seedDBIndexCache_->blockLines.size())) {
        std::ostringstream oss;
        oss << "Invalid blockId (SeedDBReader). blockId = " << blockId
            << ", blocks.size() = " << seedDBIndexCache_->blockLines.size();
        throw std::runtime_error(oss.str());
    }

    // Sequences in the block might not be stored contiguously in the file,
    // for example if a user has permuted or filtered the DB.
    // We will collect all contiguous stretches of bytes here, and then
    // fetch those parts later.
    std::vector<PacBio::Pancake::ContiguousFilePart> contiguousParts =
        GetContiguousParts(seedDBIndexCache_, blockId);

    int64_t totalBytes = 0;
    for (const auto& part : contiguousParts) {
        totalBytes += (part.endOffset - part.startOffset);
    }
    const int64_t totalItems = totalBytes / 16;

    // Preallocate the space for the loaded data.
    std::vector<SeedDB::SeedRaw> ret(totalItems);

    OpenFileHandler fh;
    int64_t readItems = 0;
    for (const auto& part : contiguousParts) {
        // Open a new file if required.
        if (part.fileId != fh.fileId) {
            if (part.fileId < 0 ||
                part.fileId >= static_cast<int32_t>(seedDBIndexCache_->fileLines.size()))
                throw std::runtime_error("(SeedDBReaderRawBlock) Invalid fileId value: " +
                                         std::to_string(part.fileId));
            const auto& fl = seedDBIndexCache_->fileLines[part.fileId];
            const std::string actualPath =
                JoinPath(seedDBIndexCache_->indexParentFolder, fl.filename);
            fh.fp = PacBio::Pancake::OpenFile(actualPath.c_str(), "rb");
            fh.fileId = part.fileId;
            fh.pos = 0;
        }
        // Jump to a different offset in the file if required.
        if (part.startOffset != fh.pos) {
            int32_t rv = fseek(fh.fp.get(), part.startOffset, SEEK_SET);
            if (rv)
                throw std::runtime_error("(SeedDBReaderRawBlock) Could not fseek to position: " +
                                         std::to_string(part.startOffset));
            fh.pos = part.startOffset;
        }

        // Load the bytes.
        const int64_t itemsToRead = (part.endOffset - part.startOffset) / 16;
        const int64_t n = fread(&ret[readItems], sizeof(SeedDB::SeedRaw), itemsToRead, fh.fp.get());
        readItems += n;

        // Sanity check.
        if (n != itemsToRead) {
            std::ostringstream oss;
            oss << "(SeedDBReaderRawBlock) Could not read seeds for the following part: "
                << "fileId = " << part.fileId << ", offsetStart = " << part.startOffset
                << ", offsetEnd = " << part.endOffset;
            throw std::runtime_error(oss.str());
        }
    }

    return ret;
}

std::vector<ContiguousFilePart> GetContiguousParts(
    const std::shared_ptr<PacBio::Pancake::SeedDBIndexCache>& seedDBIndexCache, int32_t blockId)
{
    const auto& block = seedDBIndexCache->GetBlockLine(blockId);

    if (block.startSeqId < 0 || block.endSeqId < 0 || block.endSeqId < block.startSeqId ||
        block.startSeqId >= static_cast<int32_t>(seedDBIndexCache->seedLines.size()) ||
        block.endSeqId > static_cast<int32_t>(seedDBIndexCache->seedLines.size())) {
        std::ostringstream oss;
        oss << "The SeedDB index cache is corrupt. The block's startSeqId or endSeqId "
            << "are not valid in SeedDBIndexCache. "
            << "blockId = " << blockId << ", startSeqId = " << block.startSeqId
            << ", endSeqId = " << block.endSeqId;
        throw std::runtime_error(oss.str());
    }

    // Sequences in the block might not be stored contiguously in the file,
    // for example if a user has permuted or filtered the DB.
    // We will collect all contiguous stretches of bytes here, and then
    // fetch those parts later.
    std::vector<ContiguousFilePart> contiguousParts;

    auto AddContiguousPart = [&](const SeedDBSeedsLine& sl) {
        contiguousParts.emplace_back(
            ContiguousFilePart{sl.fileId, sl.fileOffset, sl.fileOffset + sl.numBytes});
    };

    for (int32_t ordId = block.startSeqId; ordId < block.endSeqId; ++ordId) {
        const auto& sl = seedDBIndexCache->seedLines[ordId];

        if (contiguousParts.empty()) {
            AddContiguousPart(sl);

        } else if (sl.fileId != contiguousParts.back().fileId) {
            AddContiguousPart(sl);

        } else if (sl.fileOffset == contiguousParts.back().endOffset) {
            contiguousParts.back().endOffset += sl.numBytes;

        } else if (sl.fileOffset > contiguousParts.back().endOffset ||
                   (sl.fileOffset + sl.numBytes) <= contiguousParts.back().startOffset) {
            // Allow out of order byte spans, as long as there is no overlap.
            AddContiguousPart(sl);

        } else {
            // An overlap occurred.
            const auto& last = contiguousParts.back();
            std::ostringstream oss;
            oss << "Invalid SeedsLine object in the block, overlapping other SeedsLine objects in "
                   "terms of the file offset. Last ContiguousFilePart span: {"
                << last.fileId << ", " << last.startOffset << ", " << last.endOffset
                << "}, last SeedsLine: {" << sl.seqId << ", " << sl.header << ", " << sl.fileId
                << ", " << sl.fileOffset << ", " << sl.numBytes << ", " << sl.numBases << ", "
                << sl.numSeeds << "}.";
            throw std::runtime_error(oss.str());
        }
    }

    return contiguousParts;
}

}  // namespace Pancake
}  // namespace PacBio
