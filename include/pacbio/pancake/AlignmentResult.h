// Author: Ivan Sovic

#ifndef PANCAKE_ALIGNMENT_RESULT_H
#define PANCAKE_ALIGNMENT_RESULT_H

#include <pbbam/Cigar.h>

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

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_ALIGNMENT_RESULT_H
