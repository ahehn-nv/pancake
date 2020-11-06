// Author: Ivan Sovic

#ifndef PANCAKE_ALIGNER_BASE_H
#define PANCAKE_ALIGNER_BASE_H

#include <pacbio/pancake/AlignmentParameters.h>
#include <pacbio/pancake/AlignmentResult.h>
#include <cstdint>
#include <memory>
#include <vector>

namespace PacBio {
namespace Pancake {

class AlignerBase
{
public:
    virtual ~AlignerBase() = default;
    virtual AlignmentResult Global(const char* qseq, int64_t qlen, const char* tseq,
                                   int64_t tlen) = 0;
    virtual AlignmentResult Extend(const char* qseq, int64_t qlen, const char* tseq,
                                   int64_t tlen) = 0;
};

typedef std::shared_ptr<AlignerBase> AlignerBasePtr;

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_MAPPER_CLR_H
