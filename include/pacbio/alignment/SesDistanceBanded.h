// Author: Ivan Sovic

#ifndef PANCAKE_ALIGNMENT_SES_DISTANCE_BANDED_H
#define PANCAKE_ALIGNMENT_SES_DISTANCE_BANDED_H

#include <cstdint>
#include <limits>
#include <sstream>

namespace PacBio {
namespace Pancake {
namespace Alignment {

const int32_t MINUS_INF =
    std::numeric_limits<int32_t>::min() + 10000;  // We need a margin to avoid overflows.

class SesResults
{
public:
    int32_t lastQueryPos = 0;
    int32_t lastTargetPos = 0;
    int32_t diffs = 0;
    int32_t minK = 0;
    int32_t maxK = 0;
    bool valid = false;

    bool operator==(const SesResults& b) const
    {
        return lastQueryPos == b.lastQueryPos && lastTargetPos == b.lastTargetPos &&
               diffs == b.diffs && minK == b.minK && maxK == b.maxK && valid == b.valid;
    }
    friend std::ostream& operator<<(std::ostream& os, const SesResults& r);
};
inline std::ostream& operator<<(std::ostream& os, const SesResults& a)
{
    os << "lastQueryPos = " << a.lastQueryPos << ", lastTargetPos = " << a.lastTargetPos
       << ", diffs = " << a.diffs << ", minK = " << a.minK << ", maxK = " << a.maxK
       << ", valid = " << a.valid;
    return os;
}

SesResults SESDistanceBanded(const char* query, size_t queryLen, const char* target,
                             size_t targetLen, int32_t maxDiffs, int32_t bandwidth);
}
}
}

#endif  // PANCAKE_ALIGNMENT_SES_DISTANCE_BANDED_H
