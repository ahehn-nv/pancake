// Author: Ivan Sovic

#ifndef PANCAKE_ALIGNMENT_SEEDED_H
#define PANCAKE_ALIGNMENT_SEEDED_H

#include <pacbio/pancake/AlignerFactory.h>
#include <pacbio/pancake/DPChain.h>
#include <pacbio/pancake/Overlap.h>
#include <pacbio/pancake/Seed.h>
#include <cstdint>
#include <memory>
#include <vector>

namespace PacBio {
namespace Pancake {

enum class RegionType
{
    FRONT,
    BACK,
    GLOBAL,
};

inline std::string RegionTypeToString(const RegionType& type)
{
    if (type == RegionType::FRONT) {
        return "FRONT";
    } else if (type == RegionType::BACK) {
        return "BACK";
    } else if (type == RegionType::GLOBAL) {
        return "GLOBAL";
    }
    return "UNKNOWN";
}

class AlignmentRegion
{
public:
    int32_t qStart = 0;
    int32_t qSpan = 0;
    int32_t tStart = 0;
    int32_t tSpan = 0;
    bool queryRev = false;
    RegionType type = RegionType::GLOBAL;
    int32_t regionId = -1;

    bool operator==(const AlignmentRegion& b) const
    {
        return qStart == b.qStart && qSpan == b.qSpan && tStart == b.tStart && tSpan == b.tSpan &&
               queryRev == b.queryRev && type == b.type && regionId == b.regionId;
    }
};
inline std::ostream& operator<<(std::ostream& os, const AlignmentRegion& b)
{
    os << "qStart = " << b.qStart << ", qSpan = " << b.qSpan << ", tStart = " << b.tStart
       << ", tSpan = " << b.tSpan << ", queryRev = " << (b.queryRev ? "true" : "false")
       << ", type = " << RegionTypeToString(b.type) << ", regionId = " << b.regionId;
    return os;
}

class AlignRegionsGenericResult
{
public:
    Data::Cigar cigar;
    int32_t offsetFrontQuery = 0;
    int32_t offsetBackQuery = 0;
    int32_t offsetFrontTarget = 0;
    int32_t offsetBackTarget = 0;
};

std::vector<AlignmentRegion> ExtractAlignmentRegions(const std::vector<SeedHit>& inSortedHits,
                                                     int32_t qLen, int32_t tLen, bool isRev,
                                                     int32_t minAlignmentSpan,
                                                     int32_t maxFlankExtensionDist,
                                                     double flankExtensionFactor);

AlignmentResult AlignSingleRegion(const char* targetSeq, int32_t targetLen, const char* querySeqFwd,
                                  const char* querySeqRev, int32_t queryLen,
                                  AlignerBasePtr& alignerGlobal, AlignerBasePtr& alignerExt,
                                  const AlignmentRegion& region);

AlignRegionsGenericResult AlignRegionsGeneric(const char* targetSeq, const int32_t targetLen,
                                              const char* queryFwd, const char* queryRev,
                                              const int32_t queryLen,
                                              const std::vector<AlignmentRegion>& regions,
                                              AlignerBasePtr& alignerGlobal,
                                              AlignerBasePtr& alignerExt);

OverlapPtr AlignmentSeeded(const OverlapPtr& ovl, const std::vector<SeedHit>& sortedHits,
                           const char* targetSeq, const int32_t targetLen, const char* queryFwd,
                           const char* queryRev, const int32_t queryLen, int32_t minAlignmentSpan,
                           int32_t maxFlankExtensionDist, double flankExtensionFactor,
                           AlignerBasePtr& alignerGlobal, AlignerBasePtr& alignerExt);

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_ALIGNMENT_SEEDED_H
