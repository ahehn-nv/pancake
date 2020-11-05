// Authors: Ivan Sovic

#include <pacbio/pancake/SeedHit.h>
#include <sstream>

namespace PacBio {
namespace Pancake {

void CalcHitCoverage(const std::vector<SeedHit>& hits, int32_t seedLen, int32_t hitsBegin,
                     int32_t hitsEnd, int32_t& coveredBasesQuery, int32_t& coveredBasesTarget)
{
    /*
      Expects the seed hits to be sorted!
    */
    coveredBasesQuery = coveredBasesTarget = 0;

    if (hits.size() == 0 || hitsBegin >= static_cast<int32_t>(hits.size()) ||
        (hitsEnd - hitsBegin) <= 0) {
        return;
    }

    // Add the left part of the tile to the covered bases.
    hitsEnd = std::min(hitsEnd, (int32_t)hits.size());
    for (int32_t i = (hitsBegin + 1); i < hitsEnd; i++) {
        if (hits[i].queryPos < hits[i - 1].queryPos) {
            std::ostringstream oss;
            oss << "Invalid seed hit ordering, hits are either not sorted properly or not "
                   "monotonically increasing in both query and target coordinates. "
                << "hits[i - 1].queryPos = " << hits[i - 1].queryPos
                << ", hits[i - 1].targetPos = " << hits[i - 1].targetPos
                << ", hits[i].queryPos = " << hits[i].queryPos
                << ", hits[i].targetPos = " << hits[i].targetPos << ", i = " << i;
            oss << "\n";
            throw std::runtime_error(oss.str());
        }
        coveredBasesQuery += std::min(seedLen, (hits[i].queryPos - hits[i - 1].queryPos));
        coveredBasesTarget += std::min(seedLen, (hits[i].targetPos - hits[i - 1].targetPos));
    }

    // The last seed needs to be covered fully.
    coveredBasesQuery += seedLen;
    coveredBasesTarget += seedLen;
}

}  // namespace Pancake
}  // namespace PacBio
