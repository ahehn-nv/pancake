// Author: Ivan Sovic

#ifndef PANCAKE_CONTIGUOUS_FILE_PART_H
#define PANCAKE_CONTIGUOUS_FILE_PART_H

#include <cstdint>

namespace PacBio {
namespace Pancake {

class ContiguousFilePart
{
public:
    int32_t fileId = 0;
    int64_t startOffset = 0;
    int64_t endOffset = 0;
    int32_t startId = 0;
    int32_t endId = 0;

    bool operator==(const ContiguousFilePart& b) const
    {
        return fileId == b.fileId && startOffset == b.startOffset && endOffset == b.endOffset &&
               startId == b.startId && endId == b.endId;
    }
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_CONTIGUOUS_FILE_PART_H