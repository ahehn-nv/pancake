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

class AlignmentRegion
{
public:
    const char* qSeq = NULL;
    const char* tSeq = NULL;
    int32_t qStart = 0;
    int32_t qSpan = 0;
    int32_t tStart = 0;
    int32_t tSpan = 0;
    bool semiglobal = false;
};
class AlignmentRegionWithSeq
{
public:
    std::string qSeq;
    std::string tSeq;
};

class RegionsToAlign
{
public:
    AlignmentRegionWithSeq frontSemiglobal;
    std::vector<AlignmentRegion> internalGlobal;
    AlignmentRegion backSemiglobal;
    int32_t actualQueryStart = 0;
    int32_t actualQueryEnd = 0;
    int32_t actualTargetStart = 0;
    int32_t actualTargetEnd = 0;
};

class RegionsToAlignResults
{
public:
    std::vector<Data::Cigar> cigarChunks;

    // The following 4 values are the offsets produced by semiglobal alignment
    // from the actualQueryStart/actualQueryEnd/actualTargetStart/actualTargetEnd coordinates.
    int32_t offsetFrontQuery = 0;
    int32_t offsetFrontTarget = 0;
    int32_t offsetBackQuery = 0;
    int32_t offsetBackTarget = 0;
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
