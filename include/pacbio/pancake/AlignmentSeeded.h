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
    int32_t regionId = 0;

    bool operator==(const AlignmentRegion& b) const
    {
        return qStart == b.qStart && qSpan == b.qSpan && tStart == b.tStart && tSpan == b.tSpan &&
               queryRev == b.queryRev && type == b.type;
    }
};

class RegionsToAlign
{
public:
    std::vector<AlignmentRegion> regions;
    int32_t actualQueryStart = 0;
    int32_t actualQueryEnd = 0;
    int32_t actualTargetStart = 0;
    int32_t actualTargetEnd = 0;

    const char* targetSeq = NULL;
    const char* querySeqFwd = NULL;
    const char* querySeqRev = NULL;
    int32_t targetLen = 0;
    int32_t queryLen = 0;

    bool operator==(const RegionsToAlign& b) const
    {
        return regions == b.regions && actualQueryStart == b.actualQueryStart &&
               actualQueryEnd == b.actualQueryEnd && actualTargetStart == b.actualTargetStart &&
               actualTargetEnd == b.actualTargetEnd && targetLen == b.targetLen &&
               queryLen == b.queryLen;
    }
};

class AlignedRegion
{
public:
    Data::Cigar cigar;
    int32_t qSpan = 0;
    int32_t tSpan = 0;

    bool operator==(const AlignedRegion& b) const
    {
        return cigar == b.cigar && qSpan == b.qSpan && tSpan == b.tSpan;
    }
};

class RegionsToAlignResults
{
public:
    std::vector<AlignedRegion> alignedRegions;

    // The following 4 values are the offsets produced by semiglobal alignment
    // from the actualQueryStart/actualQueryEnd/actualTargetStart/actualTargetEnd coordinates.
    int32_t offsetFrontQuery = 0;
    int32_t offsetFrontTarget = 0;
    int32_t offsetBackQuery = 0;
    int32_t offsetBackTarget = 0;

    bool operator==(const RegionsToAlignResults& b) const
    {
        return alignedRegions == b.alignedRegions && offsetFrontQuery == b.offsetFrontQuery &&
               offsetFrontTarget == b.offsetFrontTarget && offsetBackQuery == b.offsetBackQuery &&
               offsetBackTarget == b.offsetBackTarget;
    }
};

RegionsToAlign ExtractAlignmentRegions(const std::vector<SeedHit>& inSortedHits,
                                       const char* querySeqFwd, const char* querySeqRev,
                                       int32_t qLen, const char* targetSeq, int32_t tLen,
                                       bool isRev, int32_t minAlignmentSpan,
                                       int32_t maxFlankExtensionDist);

RegionsToAlignResults AlignRegionsGeneric(const RegionsToAlign& regions,
                                          AlignerBasePtr& alignerGlobal,
                                          AlignerBasePtr& alignerExt);

OverlapPtr AlignmentSeeded(const OverlapPtr& ovl, const std::vector<SeedHit>& sortedHits,
                           const char* targetSeq, const int32_t targetLen, const char* queryFwd,
                           const char* queryRev, const int32_t queryLen, int32_t minAlignmentSpan,
                           int32_t maxFlankExtensionDist, AlignerBasePtr& alignerGlobal,
                           AlignerBasePtr& alignerExt);

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_ALIGNMENT_SEEDED_H
