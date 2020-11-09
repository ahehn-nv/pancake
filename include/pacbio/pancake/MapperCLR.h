// Author: Ivan Sovic

#ifndef PANCAKE_MAPPER_CLR_H
#define PANCAKE_MAPPER_CLR_H

#include <lib/intervaltree/IntervalTree.h>
#include <pacbio/overlaphifi/OverlapHifiSettings.h>
#include <pacbio/pancake/AlignerBase.h>
#include <pacbio/pancake/AlignerFactory.h>
#include <pacbio/pancake/AlignerKSW2.h>
#include <pacbio/pancake/AlignerWFA.h>
#include <pacbio/pancake/AlignmentParameters.h>
#include <pacbio/pancake/DPChain.h>
#include <pacbio/pancake/FastaSequenceCached.h>
#include <pacbio/pancake/Overlap.h>
#include <pacbio/pancake/Range.h>
#include <pacbio/pancake/Seed.h>
#include <pacbio/pancake/SeedDBIndexCache.h>
#include <pacbio/pancake/SeedIndex.h>
#include <pacbio/pancake/SeqDBReaderCachedBlock.h>
#include <pacbio/pancake/SequenceSeedsCached.h>
#include <pacbio/util/CommonTypes.h>
#include <cstdint>
#include <memory>
#include <unordered_map>
#include <vector>

namespace PacBio {
namespace Pancake {

using IntervalTreeInt32 =
    interval_tree::IntervalTree<int32_t,
                                int32_t>;  // First: interval scalar type, Second: value type.
using IntervalVectorInt32 = IntervalTreeInt32::interval_vector;
using IntervalInt32 = IntervalTreeInt32::interval;

// clang-format off
class MapperCLRSettings
{
public:
    // Mapping.
    int32_t chainMaxSkip = 25;
    int32_t chainMaxPredecessors = 500;
    int32_t chainBandwidth = 500;
    int32_t maxGap = 10000;
    int32_t minNumSeeds = 3;
    int32_t minCoveredBases = 0;
    int32_t minDPScore = 40;
    double secondaryAllowedOverlapFractionQuery = 0.50;
    double secondaryAllowedOverlapFractionTarget = 0.50;
    double secondaryMinScoreFraction = 0.80;
    bool useLIS = true;

    // Indexing.
    PacBio::Pancake::SeedDB::SeedDBParameters seedParams;
    PacBio::Pancake::SeedDB::SeedDBParameters seedParamsFallback;
    double freqPercentile = 0.000;

    // Alignment.
    bool align = true;
    int32_t maxFlankExtensionDist = 10000;      // Maximum length of the query/target sequences to consider when aligning flanks.
    int32_t minAlignmentSpan = 200;             // If two seeds are closer than this, take the next seed. (Unless there are only 2 seeds left.)
    AlignerType alignerTypeGlobal = AlignerType::KSW2;
    AlignmentParameters alnParamsGlobal;
    AlignerType alignerTypeExt = AlignerType::KSW2;
    AlignmentParameters alnParamsExt;

    // Other.
    bool skipSymmetricOverlaps = false;
    int32_t minQueryLen = 50;
    int32_t bestNSecondary = 0;

    // bool oneHitPerTarget = false;
    // int32_t maxSeedDistance = 5000;
    // int32_t minMappedLength = 1000;
    // double minIdentity = 75.0;
};
// clang-format on

inline std::ostream& operator<<(std::ostream& out, const MapperCLRSettings& a)
{
    out << "chainMaxSkip = " << a.chainMaxSkip << "\n"
        << "chainMaxPredecessors = " << a.chainMaxPredecessors << "\n"
        << "chainBandwidth = " << a.chainBandwidth << "\n"
        << "maxGap = " << a.maxGap << "\n"
        << "minNumSeeds = " << a.minNumSeeds << "\n"
        << "minCoveredBases = " << a.minCoveredBases << "\n"
        << "minDPScore = " << a.minDPScore << "\n"
        << "useLIS = " << a.useLIS << "\n"

        << "secondaryAllowedOverlapFractionQuery = " << a.secondaryAllowedOverlapFractionQuery
        << "\n"
        << "secondaryAllowedOverlapFractionTarget = " << a.secondaryAllowedOverlapFractionTarget
        << "\n"
        << "secondaryMinScoreFraction = " << a.secondaryMinScoreFraction << "\n"

        << "align = " << a.align << "\n"
        << "maxFlankExtensionDist = " << a.maxFlankExtensionDist << "\n"
        << "minAlignmentSpan = " << a.minAlignmentSpan << "\n"
        << "alignerTypeGlobal = " << AlignerTypeToString(a.alignerTypeGlobal) << "\n"
        << "alnParamsGlobal:\n"
        << a.alnParamsGlobal << "alignerTypeExt = " << AlignerTypeToString(a.alignerTypeExt) << "\n"
        << "alnParamsExt:\n"
        << a.alnParamsExt

        << "seedParams.KmerSize = " << a.seedParams.KmerSize << "\n"
        << "seedParams.MinimizerWindow = " << a.seedParams.MinimizerWindow << "\n"
        << "seedParams.Spacing = " << a.seedParams.Spacing << "\n"
        << "seedParams.UseHPC = " << a.seedParams.UseHPC << "\n"
        << "seedParams.UseHPCForSeedsOnly = " << a.seedParams.UseHPCForSeedsOnly << "\n"
        << "seedParams.MaxHPCLen = " << a.seedParams.MaxHPCLen << "\n"
        << "seedParams.UseRC = " << a.seedParams.UseRC << "\n"

        << "seedParamsFallback.KmerSize = " << a.seedParamsFallback.KmerSize << "\n"
        << "seedParamsFallback.MinimizerWindow = " << a.seedParamsFallback.MinimizerWindow << "\n"
        << "seedParamsFallback.Spacing = " << a.seedParamsFallback.Spacing << "\n"
        << "seedParamsFallback.UseHPC = " << a.seedParamsFallback.UseHPC << "\n"
        << "seedParamsFallback.UseHPCForSeedsOnly = " << a.seedParamsFallback.UseHPCForSeedsOnly
        << "\n"
        << "seedParamsFallback.MaxHPCLen = " << a.seedParamsFallback.MaxHPCLen << "\n"
        << "seedParamsFallback.UseRC = " << a.seedParamsFallback.UseRC << "\n"

        << "freqPercentile = " << a.freqPercentile << "\n"

        << "skipSymmetricOverlaps = " << a.skipSymmetricOverlaps << "\n"
        << "minQueryLen = " << a.minQueryLen << "\n"
        << "bestNSecondary = " << a.bestNSecondary << "\n";

    return out;
}

struct ChainedRegion
{
    ChainedHits chain;
    OverlapPtr mapping;
    int32_t priority = 0;  // Priority 0 means primary alignment.
    bool isSupplementary = false;
};

class MapperCLRResult
{
public:
    std::vector<std::unique_ptr<ChainedRegion>> mappings;
};

class MapperCLR
{
public:
    MapperCLR(const MapperCLRSettings& settings);
    ~MapperCLR();

    std::vector<MapperCLRResult> Map(const std::vector<std::string>& targetSeqs,
                                     const std::vector<std::string>& querySeqs);

    MapperCLRResult Map(const std::vector<std::string>& targetSeqs,
                        const PacBio::Pancake::SeedIndex& index, const std::string& querySeq,
                        const std::vector<PacBio::Pancake::Int128t>& querySeeds,
                        const int32_t queryId, int64_t freqCutoff);

private:
    MapperCLRSettings settings_;
    AlignerBasePtr alignerGlobal_;
    AlignerBasePtr alignerExt_;

    static OverlapPtr MakeOverlap_(const std::vector<SeedHit>& sortedHits, int32_t queryId,
                                   int32_t queryLen, const PacBio::Pancake::SeedIndex& index,
                                   int32_t beginId, int32_t endId, int32_t minTargetPosId,
                                   int32_t maxTargetPosId);

    static std::vector<Range> DiagonalGroup_(const std::vector<SeedHit>& sortedHits,
                                             int32_t chainBandwidth, bool overlappingWindows);

    static std::vector<std::unique_ptr<ChainedRegion>> ChainAndMakeOverlap_(
        const PacBio::Pancake::SeedIndex& index, const std::vector<SeedHit>& hits,
        const std::vector<PacBio::Pancake::Range>& hitGroups, int32_t queryId, int32_t queryLen,
        int32_t chainMaxSkip, int32_t chainMaxPredecessors, int32_t maxGap, int32_t chainBandwidth,
        int32_t minNumSeeds, int32_t minCoveredBases, int32_t minDPScore, bool useLIS);

    static std::vector<std::unique_ptr<ChainedRegion>> ReChainSeedHits_(
        const std::vector<std::unique_ptr<ChainedRegion>>& chainedRegions,
        const PacBio::Pancake::SeedIndex& index, int32_t queryId, int32_t queryLen,
        int32_t chainMaxSkip, int32_t chainMaxPredecessors, int32_t maxGap, int32_t chainBandwidth,
        int32_t minNumSeeds, int32_t minCoveredBases, int32_t minDPScore);

    static void WriteSeedHits_(const std::string& outPath, const std::vector<SeedHit>& hits,
                               size_t hitsStart, size_t hitsEnd, int32_t hitsId,
                               const std::string& queryName, int64_t queryLength,
                               const std::string& targetName, int64_t targetLength, bool append);

    static std::vector<Range> GroupByTargetAndStrand_(const std::vector<SeedHit>& sortedHits);

    static void LongMergeChains_(std::vector<std::unique_ptr<ChainedRegion>>& chainedRegions,
                                 int32_t maxGap);

    static void WrapFlagSecondaryAndSupplementary_(
        std::vector<std::unique_ptr<ChainedRegion>>& allChainedRegions,
        double secondaryAllowedOverlapFractionQuery, double secondaryAllowedOverlapFractionTarget,
        double secondaryMinScoreFraction);
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_MAPPER_CLR_H
