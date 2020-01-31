// Author: Ivan Sovic

#ifndef PANCAKE_SEEDDB_READER_RAW_BLOCK_H
#define PANCAKE_SEEDDB_READER_RAW_BLOCK_H

#include <pacbio/seeddb/Seed.h>
#include <pacbio/seeddb/SeedDBIndexCache.h>
#include <pacbio/seqdb/Util.h>
#include <memory>
#include <string>
#include <vector>

namespace PacBio {
namespace Pancake {

class SeedDBReaderRawBlock
{
public:
    SeedDBReaderRawBlock(const std::shared_ptr<PacBio::Pancake::SeedDBIndexCache>& seedDBCache);
    ~SeedDBReaderRawBlock() = default;

    std::vector<SeedDB::SeedRaw> GetBlock(int32_t blockId) const;

private:
    using FilePtr = std::unique_ptr<FILE, FileDeleter>;
    class OpenFileHandler
    {
    public:
        FilePtr fp = nullptr;
        int32_t fileId = -1;
        int64_t pos = -1;
    };

    const std::shared_ptr<PacBio::Pancake::SeedDBIndexCache> seedDBIndexCache_;
};

class ContiguousFilePart
{
public:
    int32_t fileId = 0;
    int64_t startOffset = 0;
    int64_t endOffset = 0;

    bool operator==(const ContiguousFilePart& b) const
    {
        return fileId == b.fileId && startOffset == b.startOffset && endOffset == b.endOffset;
    }
};
std::vector<ContiguousFilePart> GetContiguousParts(
    const std::shared_ptr<PacBio::Pancake::SeedDBIndexCache>& seedDBIndexCache, int32_t blockId);

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_SEEDDB_READER_RAW_BLOCK_H