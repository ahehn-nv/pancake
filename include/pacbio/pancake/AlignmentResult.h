// Author: Ivan Sovic

#ifndef PANCAKE_ALIGNMENT_RESULT_H
#define PANCAKE_ALIGNMENT_RESULT_H

#include <pbbam/Cigar.h>
#include <ostream>

namespace PacBio {
namespace Pancake {

class AlignmentResult
{
public:
    PacBio::BAM::Cigar cigar;
    int32_t lastQueryPos = 0;
    int32_t lastTargetPos = 0;
    int32_t maxQueryPos = 0;
    int32_t maxTargetPos = 0;
    bool valid = false;
    int32_t score = 0;
    int32_t maxScore = 0;
    bool zdropped = false;
};

inline std::ostream& operator<<(std::ostream& os, const AlignmentResult& b)
{
    os << "valid = " << (b.valid ? "true" : "false") << ", score = " << b.score
       << ", maxScore = " << b.maxScore << ", zdropped = " << b.zdropped
       << ", lastQueryPos = " << b.lastQueryPos << ", lastTargetPos = " << b.lastTargetPos
       << ", maxQueryPos = " << b.maxQueryPos << ", maxTargetPos = " << b.maxTargetPos
       << ", CIGAR: " << b.cigar.ToStdString();
    return os;
}

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_ALIGNMENT_RESULT_H
