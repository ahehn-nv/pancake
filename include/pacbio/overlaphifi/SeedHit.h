// Author: Ivan Sovic

#ifndef PANCAKE_OVERLAPHIFI_SEEDHIT_H
#define PANCAKE_OVERLAPHIFI_SEEDHIT_H

#include <cstdint>

namespace PacBio {
namespace Pancake {

class SeedHit
{
public:
    SeedHit() = default;
    ~SeedHit() = default;
    __int128 PackTo128() const
    {
        __int128 ret = 0;
        ret = (static_cast<__int128>(targetId) << 97) | (static_cast<__int128>(targetRev) << 96) |
              (static_cast<__int128>(targetPos) << 64) | (static_cast<__int128>(flags) << 32) |
              static_cast<__int128>(queryPos);
        return ret;
    }
    bool operator<(const SeedHit& b) const { return this->PackTo128() < b.PackTo128(); }
    bool operator==(const SeedHit& b) const
    {
        return targetId == b.targetId && targetRev == b.targetRev && targetPos == b.targetPos &&
               flags == b.flags && queryPos == b.queryPos;
    }

    int32_t targetId : 31;
    bool targetRev : 1;
    int32_t targetPos : 32;
    int32_t flags : 32;
    int32_t queryPos : 32;
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_OVERLAPHIFI_SEEDHIT_H