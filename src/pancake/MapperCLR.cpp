// Authors: Ivan Sovic

#include <lib/kxsort/kxsort.h>
#include <pacbio/alignment/AlignmentTools.h>
#include <pacbio/alignment/DiffCounts.h>
#include <pacbio/alignment/SesDistanceBanded.h>
#include <pacbio/pancake/AlignmentSeeded.h>
#include <pacbio/pancake/MapperCLR.h>
#include <pacbio/pancake/Minimizers.h>
#include <pacbio/pancake/OverlapWriterBase.h>
#include <pacbio/pancake/Secondary.h>
#include <pacbio/pancake/SeedHitWriter.h>
#include <pacbio/util/RunLengthEncoding.h>
#include <pacbio/util/TicToc.h>
#include <pacbio/util/Util.h>
#include <pbcopper/logging/Logging.h>
#include <pbcopper/third-party/edlib.h>
#include <algorithm>
#include <array>
#include <iostream>
#include <lib/istl/lis.hpp>
#include <sstream>
#include <tuple>

namespace PacBio {
namespace Pancake {

// #define PANCAKE_MAP_CLR_DEBUG
// #define PANCAKE_MAP_CLR_DEBUG_2
// #define PANCAKE_WRITE_SCATTERPLOT
// #define PANCAKE_MAP_CLR_DEBUG_ALIGN

MapperCLR::MapperCLR(const MapperCLRSettings& settings)
    : settings_{settings}, alignerGlobal_(nullptr), alignerExt_(nullptr)
{
    alignerGlobal_ = AlignerFactory(settings.alignerTypeGlobal, settings.alnParamsGlobal);
    alignerExt_ = AlignerFactory(settings.alignerTypeExt, settings.alnParamsExt);
}

MapperCLR::~MapperCLR() = default;

std::vector<MapperBaseResult> MapperCLR::MapAndAlign(const std::vector<std::string>& targetSeqs,
                                                     const std::vector<std::string>& querySeqs)
{
    std::vector<FastaSequenceCached> targetSeqsCached;
    for (int32_t i = 0; i < static_cast<int32_t>(targetSeqs.size()); ++i) {
        targetSeqsCached.emplace_back(
            FastaSequenceCached(std::to_string(i), targetSeqs[i].c_str(), targetSeqs[i].size(), i));
    }

    std::vector<FastaSequenceCached> querySeqsCached;
    for (int32_t i = 0; i < static_cast<int32_t>(querySeqs.size()); ++i) {
        querySeqsCached.emplace_back(
            FastaSequenceCached(std::to_string(i), querySeqs[i].c_str(), querySeqs[i].size(), i));
    }

    return MapAndAlign(targetSeqsCached, querySeqsCached);
}

std::vector<MapperBaseResult> MapperCLR::MapAndAlign(const std::vector<FastaSequenceId>& targetSeqs,
                                                     const std::vector<FastaSequenceId>& querySeqs)
{
    std::vector<FastaSequenceCached> targetSeqsCached;
    for (int32_t i = 0; i < static_cast<int32_t>(targetSeqs.size()); ++i) {
        targetSeqsCached.emplace_back(
            FastaSequenceCached(std::to_string(i), targetSeqs[i].Bases().c_str(),
                                targetSeqs[i].Bases().size(), targetSeqs[i].Id()));
    }

    std::vector<FastaSequenceCached> querySeqsCached;
    for (int32_t i = 0; i < static_cast<int32_t>(querySeqs.size()); ++i) {
        querySeqsCached.emplace_back(
            FastaSequenceCached(std::to_string(i), querySeqs[i].Bases().c_str(),
                                querySeqs[i].Bases().size(), querySeqs[i].Id()));
    }

    return MapAndAlign(targetSeqsCached, querySeqsCached);
}

std::vector<MapperBaseResult> MapperCLR::MapAndAlign(
    const std::vector<FastaSequenceCached>& targetSeqs,
    const std::vector<FastaSequenceCached>& querySeqs)
{
    // Construct the index.
    std::vector<PacBio::Pancake::Int128t> seeds;
    std::vector<int32_t> sequenceLengths;
    const auto& seedParams = settings_.seedParams;
    SeedDB::GenerateMinimizers(seeds, sequenceLengths, targetSeqs, seedParams.KmerSize,
                               seedParams.MinimizerWindow, seedParams.Spacing, seedParams.UseRC,
                               seedParams.UseHPCForSeedsOnly, seedParams.MaxHPCLen);
    std::unique_ptr<SeedIndex> seedIndex =
        std::make_unique<SeedIndex>(settings_.seedParams, sequenceLengths, std::move(seeds));

    // Calculate the seed frequency statistics, needed for the cutoff.
    int64_t freqMax = 0;
    int64_t freqCutoff = 0;
    double freqAvg = 0.0;
    double freqMedian = 0.0;
    seedIndex->ComputeFrequencyStats(settings_.freqPercentile, freqMax, freqAvg, freqMedian,
                                     freqCutoff);

    // Construct the fallback index only if required
    int64_t freqCutoffFallback = 0;
    std::unique_ptr<SeedIndex> seedIndexFallback = nullptr;
    if (settings_.seedParamsFallback != settings_.seedParams) {
        std::vector<PacBio::Pancake::Int128t> seedsFallback;
        const auto& seedParamsFallback = settings_.seedParamsFallback;
        SeedDB::GenerateMinimizers(seedsFallback, sequenceLengths, targetSeqs,
                                   seedParamsFallback.KmerSize, seedParamsFallback.MinimizerWindow,
                                   seedParamsFallback.Spacing, seedParamsFallback.UseRC,
                                   seedParamsFallback.UseHPCForSeedsOnly,
                                   seedParamsFallback.MaxHPCLen);
        seedIndexFallback = std::make_unique<SeedIndex>(settings_.seedParamsFallback,
                                                        sequenceLengths, std::move(seedsFallback));
        seedIndexFallback->ComputeFrequencyStats(settings_.freqPercentile, freqMax, freqAvg,
                                                 freqMedian, freqCutoffFallback);
    }

    // Run mapping for each query.
    std::vector<MapperBaseResult> results;
    for (int32_t i = 0; i < static_cast<int32_t>(querySeqs.size()); ++i) {
        const auto& query = querySeqs[i];
        int32_t queryId = query.Id();

#ifdef PANCAKE_MAP_CLR_DEBUG
        PBLOG_INFO << "[i = " << i << "] New query: queryId = " << queryId
                   << ", length = " << query.size();
#endif

        std::vector<PacBio::Pancake::Int128t> querySeeds;
        int32_t seqLen = query.size();
        const uint8_t* seq = reinterpret_cast<const uint8_t*>(query.data());
        int rv = SeedDB::GenerateMinimizers(
            querySeeds, seq, seqLen, 0, queryId, settings_.seedParams.KmerSize,
            settings_.seedParams.MinimizerWindow, settings_.seedParams.Spacing,
            settings_.seedParams.UseRC, settings_.seedParams.UseHPCForSeedsOnly,
            settings_.seedParams.MaxHPCLen);
        if (rv)
            throw std::runtime_error("Generating minimizers failed for the query sequence i = " +
                                     std::to_string(i) + ", id = " + std::to_string(queryId));

        auto queryResults =
            MapAndAlign(targetSeqs, *seedIndex, query, querySeeds, queryId, freqCutoff);

        if (queryResults.mappings.empty() && seedIndexFallback != nullptr) {
            rv = SeedDB::GenerateMinimizers(
                querySeeds, seq, seqLen, 0, queryId, settings_.seedParamsFallback.KmerSize,
                settings_.seedParamsFallback.MinimizerWindow, settings_.seedParamsFallback.Spacing,
                settings_.seedParamsFallback.UseRC, settings_.seedParamsFallback.UseHPCForSeedsOnly,
                settings_.seedParamsFallback.MaxHPCLen);
            if (rv)
                throw std::runtime_error(
                    "Generating minimizers failed for the query sequence, id = " +
                    std::to_string(queryId));

            queryResults = MapAndAlign(targetSeqs, *seedIndexFallback, query, querySeeds, queryId,
                                       freqCutoffFallback);
        }

        for (const auto& m : queryResults.mappings) {
            m->mapping->Aid = queryId;
        }
        results.emplace_back(std::move(queryResults));

#ifdef PANCAKE_MAP_CLR_DEBUG
        PBLOG_INFO << "\n\n\n";
#endif
    }

    return results;
}

void DebugPrintChainedRegion(std::ostream& oss, int32_t regionId, const ChainedRegion& cr)
{
    oss << "[regionId " << regionId << "] chain.hits = " << cr.chain.hits.size()
        << ", chain.score = " << cr.chain.score << ", chain.covQ = " << cr.chain.coveredBasesQuery
        << ", chain.covT = " << cr.chain.coveredBasesTarget << ", priority = " << cr.priority
        << ", isSuppl = " << (cr.isSupplementary ? "true" : "false")
        << ", ovl: " << OverlapWriterBase::PrintOverlapAsM4(cr.mapping, "", "", true, false)
        << ", diagStart = " << (cr.mapping->Astart - cr.mapping->Bstart)
        << ", diagEnd = " << (cr.mapping->Aend - cr.mapping->Bend);

    // for (size_t j = 0; j < cr.chain.hits.size(); ++j) {
    //     const auto& hit = cr.chain.hits[j];
    //     std::cerr << "    [hit " << j << "] tid = " << hit.targetId
    //                 << ", trev = " << hit.targetRev << ", tpos = " << hit.targetPos
    //                 << ", qpos = " << hit.queryPos << "\n";
    // }
}

void DebugWriteChainedRegion(
#ifdef PANCAKE_MAP_CLR_DEBUG_2
    const std::vector<std::unique_ptr<ChainedRegion>>& allChainedRegions,
    const std::string& descriptor, int32_t queryId, int32_t queryLen
#else
    const std::vector<std::unique_ptr<ChainedRegion>>& /*allChainedRegions*/,
    const std::string& /*descriptor*/, int32_t /*queryId*/, int32_t /*queryLen*/
#endif
    )
{
#ifdef PANCAKE_MAP_CLR_DEBUG_2
    std::cerr << "(DebugWriteChainedRegion) " << descriptor
              << ": allChainedRegions.size() = " << allChainedRegions.size() << "\n";

    // Write seed hits after the first chaining stage.
    for (size_t i = 0; i < allChainedRegions.size(); ++i) {
        auto& region = allChainedRegions[i];
        if (region->priority > 1) {
            continue;
        }
        std::cerr << "    ";
        DebugPrintChainedRegion(std::cerr, i, *region);
        std::cerr << "\n\n";
        WriteSeedHits("temp-debug/hits-q" + std::to_string(queryId) + "-" + descriptor + ".csv",
                      region->chain.hits, 0, region->chain.hits.size(), i,
                      "query" + std::to_string(queryId), queryLen,
                      "target" + std::to_string(region->mapping->Bid), region->mapping->Blen,
                      (i > 0));
    }
#endif
}

MapperBaseResult MapperCLR::MapAndAlign(const std::vector<FastaSequenceCached>& targetSeqs,
                                        const PacBio::Pancake::SeedIndex& index,
                                        const FastaSequenceCached& querySeq,
                                        const std::vector<PacBio::Pancake::Int128t>& querySeeds,
                                        const int32_t queryId, int64_t freqCutoff)
{
    const int32_t queryLen = querySeq.size();

    // Map the query.
    TicToc ttMap;
    auto result = Map(index, querySeeds, queryLen, queryId, freqCutoff);
    ttMap.Stop();

    // Align if needed.
    TicToc ttAlign;
    if (settings_.align) {
        result = Align(targetSeqs, querySeq, result);
    }
    ttAlign.Stop();

    // Filter mappings.
    size_t numValid = 0;
    int32_t numSelectedSecondary = 0;
    for (size_t i = 0; i < result.mappings.size(); ++i) {
        const auto& region = result.mappings[i];
        if (region->mapping == nullptr) {
            continue;
        }
        if (region->priority > 1) {
            continue;
        }
        if (region->priority == 1 && numSelectedSecondary >= settings_.bestNSecondary) {
            continue;
        }

        // Since the secondary/supplementary labelling was repeated, we need to filter secondary
        // alignments again.
        if (region->priority == 1 && numSelectedSecondary < settings_.bestNSecondary) {
            ++numSelectedSecondary;
        }
        if (i != numValid) {
            std::swap(result.mappings[i], result.mappings[numValid]);
        }
        ++numValid;
    }
    result.mappings.resize(numValid);
    DebugWriteChainedRegion(result.mappings, "9-result-final", queryId, queryLen);

    return result;
}

MapperBaseResult MapperCLR::Map(const PacBio::Pancake::SeedIndex& index,
                                const std::vector<PacBio::Pancake::Int128t>& querySeeds,
                                const int32_t queryLen, const int32_t queryId, int64_t freqCutoff)
{
    // Skip short queries.
    if (queryLen < settings_.minQueryLen) {
        return {};
    }

    // Collect seed hits.
    TicToc ttCollectHits;
    std::vector<SeedHit> hits;
    index.CollectHits(&querySeeds[0], querySeeds.size(), queryLen, hits, freqCutoff);
    ttCollectHits.Stop();

    // Sort the seed hits.
    TicToc ttSortHits;
    std::sort(hits.begin(), hits.end(), [](const auto& a, const auto& b) {
        return PackSeedHitWithDiagonalToTuple(a) < PackSeedHitWithDiagonalToTuple(b);
    });
    ttSortHits.Stop();

    // Group seed hits by diagonal.
    TicToc ttDiagonalGroup;
    auto groups = DiagonalGroup(hits, settings_.chainBandwidth, true);
    ttDiagonalGroup.Stop();

    // Process each diagonal bin to get the chains.
    TicToc ttChain;
    std::vector<std::unique_ptr<ChainedRegion>> allChainedRegions = ChainAndMakeOverlap_(
        index, hits, groups, queryId, queryLen, settings_.chainMaxSkip,
        settings_.chainMaxPredecessors, settings_.maxGap, settings_.chainBandwidth,
        settings_.minNumSeeds, settings_.minCoveredBases, settings_.minDPScore, settings_.useLIS);
    ttChain.Stop();
    DebugWriteChainedRegion(allChainedRegions, "1-chain-and-make-overlap", queryId, queryLen);

    // Take the remaining regions, merge all seed hits, and rechain.
    // Needed because diagonal chaining was greedy and a wide window could have split
    // otherwise good chains. On the other hand, diagonal binning was needed for speed
    // in low-complexity regions.
    TicToc ttRechain;
    allChainedRegions =
        ReChainSeedHits_(allChainedRegions, index, queryId, queryLen, settings_.chainMaxSkip,
                         settings_.chainMaxPredecessors, settings_.maxGap, settings_.chainBandwidth,
                         settings_.minNumSeeds, settings_.minCoveredBases, settings_.minDPScore);
    ttRechain.Stop();
    DebugWriteChainedRegion(allChainedRegions, "2-rechain-hits", queryId, queryLen);

    // Sort all chains in descending order of the number of hits.
    std::sort(allChainedRegions.begin(), allChainedRegions.end(),
              [](const auto& a, const auto& b) { return a->chain.score > b->chain.score; });

    // Secondary/supplementary flagging.
    WrapFlagSecondaryAndSupplementary_(
        allChainedRegions, settings_.secondaryAllowedOverlapFractionQuery,
        settings_.secondaryAllowedOverlapFractionTarget, settings_.secondaryMinScoreFraction);
    DebugWriteChainedRegion(allChainedRegions, "3-wrap-flag-secondary-suppl-1", queryId, queryLen);

    // Merge long gaps.
    LongMergeChains_(allChainedRegions, settings_.maxGap);
    DebugWriteChainedRegion(allChainedRegions, "4-long-merge-chains", queryId, queryLen);

    // Again relabel, because some chains are longer now.
    WrapFlagSecondaryAndSupplementary_(
        allChainedRegions, settings_.secondaryAllowedOverlapFractionQuery,
        settings_.secondaryAllowedOverlapFractionTarget, settings_.secondaryMinScoreFraction);
    DebugWriteChainedRegion(allChainedRegions, "5-wrap-flag-secondary-suppl-2", queryId, queryLen);

    // Sort all chains by priority and then score.
    std::sort(allChainedRegions.begin(), allChainedRegions.end(),
              [](const std::unique_ptr<ChainedRegion>& a, const std::unique_ptr<ChainedRegion>& b) {
                  auto at = std::tuple<int32_t, bool, int32_t>(a->priority, a->isSupplementary,
                                                               a->chain.score);
                  auto bt = std::tuple<int32_t, bool, int32_t>(b->priority, b->isSupplementary,
                                                               b->chain.score);
                  return at < bt;
              });

    // Refine seed hits.
    for (size_t i = 0; i < allChainedRegions.size(); ++i) {
        auto& region = allChainedRegions[i];
        if (region->priority > 1) {
            continue;
        }
        ChainedHits newChain = RefineBadEnds(
            region->chain, settings_.alnParamsGlobal.alignBandwidth, settings_.minDPScore * 2);
        newChain = RefineChainedHits(newChain, 10, 40, settings_.maxGap / 2, 10);
        newChain = RefineChainedHits2(newChain, 30, settings_.maxGap / 2);

        std::swap(region->chain, newChain);

        region->mapping->Astart = region->chain.hits.front().queryPos;
        region->mapping->Bstart = region->chain.hits.front().targetPos;
        region->mapping->Aend = region->chain.hits.back().queryPos;
        region->mapping->Bend = region->chain.hits.back().targetPos;
        region->mapping->NumSeeds = region->chain.hits.size();
    }
    DebugWriteChainedRegion(allChainedRegions, "6-refining-seed-hits", queryId, queryLen);

    // Filter out the mappings.
    MapperBaseResult result;
    int32_t numSelectedSecondary = 0;
    for (size_t i = 0; i < allChainedRegions.size(); ++i) {
        const auto& region = allChainedRegions[i];
        if (region->priority > 1) {
            continue;
        }
        if (region->mapping == nullptr) {
            continue;
        }
        if (region->priority == 0) {
            result.mappings.emplace_back(std::move(allChainedRegions[i]));
        } else if (region->priority == 1 && numSelectedSecondary < settings_.bestNSecondary) {
            result.mappings.emplace_back(std::move(allChainedRegions[i]));
            ++numSelectedSecondary;
        }
    }
    DebugWriteChainedRegion(result.mappings, "7-result-mappings", queryId, queryLen);

#ifdef PANCAKE_MAP_CLR_DEBUG_2
    std::cerr << "All hits: hits.size() = " << hits.size() << "\n";
    std::cerr << "Diagonal groups: groups.size() = " << groups.size() << "\n";
    for (size_t i = 0; i < groups.size(); ++i) {
        const int32_t firstDiag = hits[groups[i].start].Diagonal();
        const int32_t lastDiag = hits[groups[i].end - 1].Diagonal();
        std::cerr << "[queryId = " << queryId << ", group " << i << "] start = " << groups[i].start
                  << ", end = " << groups[i].end << ", diagStart = " << firstDiag
                  << ", diagEnd = " << lastDiag << "\n";
    }

    // Write ALL seed hits.
    for (size_t i = 0; i < groups.size(); ++i) {
        const auto& g = groups[i];
        const int32_t targetId = hits[g.start].targetId;
        const int32_t targetLen = index.GetSequenceLength(targetId);
        WriteSeedHits(
            "temp-debug/hits-q" + std::to_string(queryId) + "-0-all-hits-diagonal-groupped.csv",
            hits, g.start, g.end, i, "query" + std::to_string(queryId), queryLen,
            "target" + std::to_string(targetId), targetLen, (i > 0));
    }
#endif

    return result;
}

MapperBaseResult MapperCLR::Align(const std::vector<FastaSequenceCached>& targetSeqs,
                                  const FastaSequenceCached& querySeq,
                                  const MapperBaseResult& mappingResult)
{
#ifdef PANCAKE_MAP_CLR_DEBUG_ALIGN
    std::cerr << "Aligning.\n";
    std::cerr << "alignerTypeGlobal = " << AlignerTypeToString(settings_.alignerTypeGlobal) << "\n";
    std::cerr << "alignerTypeExt = " << AlignerTypeToString(settings_.alignerTypeExt) << "\n";
#endif

    const int32_t queryLen = querySeq.size();
    const int32_t queryId = querySeq.Id();

    MapperBaseResult alignedResult;

    // Reverse the query sequence, needed for alignment.
    const std::string querySeqRev =
        PacBio::Pancake::ReverseComplement(querySeq.c_str(), 0, querySeq.size());

    TicToc ttAlign;
    for (size_t i = 0; i < mappingResult.mappings.size(); ++i) {
        // const auto& chain = result.mappings[i]->chain;
        const auto& chain = mappingResult.mappings[i]->chain;
        const auto& ovl = mappingResult.mappings[i]->mapping;
        const auto& tSeqFwd = targetSeqs[ovl->Bid];

        // Use a custom aligner to align.
        TicToc ttAlignmentSeeded;
        auto newOvl =
            AlignmentSeeded(ovl, chain.hits, tSeqFwd.c_str(), tSeqFwd.size(), &querySeq.c_str()[0],
                            &querySeqRev[0], queryLen, settings_.minAlignmentSpan,
                            settings_.maxFlankExtensionDist, alignerGlobal_, alignerExt_);
        ttAlignmentSeeded.Stop();

        auto newChainedRegion = std::make_unique<ChainedRegion>();
        newChainedRegion->chain = chain;
        newChainedRegion->mapping = std::move(newOvl);
        newChainedRegion->priority = mappingResult.mappings[i]->priority;
        newChainedRegion->isSupplementary = mappingResult.mappings[i]->isSupplementary;
        alignedResult.mappings.emplace_back(std::move(newChainedRegion));

#ifdef PANCAKE_MAP_CLR_DEBUG_ALIGN
        std::cerr << "[mapping i = " << i << ", before alignment] ovl: "
                  << OverlapWriterBase::PrintOverlapAsM4(ovl, "", "", true, true) << "\n";
        const auto& updatedOvl = alignedResult.mappings[i]->mapping;
        std::cerr << "[mapping i = " << i << ", after alignment] ovl: ";
        if (updatedOvl != nullptr) {
            std::cerr << OverlapWriterBase::PrintOverlapAsM4(updatedOvl, "", "", true, true)
                      << "\n";
        } else {
            std::cerr << "nullptr\n";
        }
        std::cerr << "ttAlignmentSeeded = " << ttAlignmentSeeded.GetCpuMillisecs() << " ms\n";
#endif
    }

    // Secondary/supplementary flagging.
    WrapFlagSecondaryAndSupplementary_(
        alignedResult.mappings, settings_.secondaryAllowedOverlapFractionQuery,
        settings_.secondaryAllowedOverlapFractionTarget, settings_.secondaryMinScoreFraction);

    ttAlign.Stop();

    DebugWriteChainedRegion(alignedResult.mappings, "8-result-after-align", queryId, queryLen);

    return alignedResult;
}

void MapperCLR::WrapFlagSecondaryAndSupplementary_(
    std::vector<std::unique_ptr<ChainedRegion>>& allChainedRegions,
    double secondaryAllowedOverlapFractionQuery, double secondaryAllowedOverlapFractionTarget,
    double secondaryMinScoreFraction)
{
    /*
     * Edits the objects in place.
    */
    // Copy the overlaps so we can satisfy the FlagSecondaryAndSupplementary API.
    std::vector<OverlapPtr> tmpOverlaps;
    for (size_t i = 0; i < allChainedRegions.size(); ++i) {
        if (allChainedRegions[i] == nullptr || allChainedRegions[i]->mapping == nullptr) {
            continue;
        }
        tmpOverlaps.emplace_back(createOverlap(allChainedRegions[i]->mapping));
    }
    // Flag the secondary and supplementary overlaps.
    std::vector<OverlapPriority> overlapPriorities = FlagSecondaryAndSupplementary(
        tmpOverlaps, secondaryAllowedOverlapFractionQuery, secondaryAllowedOverlapFractionTarget,
        secondaryMinScoreFraction);
    // Set the flags.
    for (size_t i = 0; i < allChainedRegions.size(); ++i) {
        if (allChainedRegions[i] == nullptr || allChainedRegions[i]->mapping == nullptr) {
            continue;
        }
        auto& cr = allChainedRegions[i];
        cr->mapping->IsSecondary = (overlapPriorities[i].priority > 0);
        cr->mapping->IsSupplementary = overlapPriorities[i].isSupplementary;
        cr->priority = overlapPriorities[i].priority;
        cr->isSupplementary = overlapPriorities[i].isSupplementary;
    }
}

std::vector<std::unique_ptr<ChainedRegion>> MapperCLR::ReChainSeedHits_(
    const std::vector<std::unique_ptr<ChainedRegion>>& chainedRegions,
    const PacBio::Pancake::SeedIndex& index, int32_t queryId, int32_t queryLen,
    int32_t chainMaxSkip, int32_t chainMaxPredecessors, int32_t maxGap, int32_t chainBandwidth,
    int32_t minNumSeeds, int32_t minCoveredBases, int32_t minDPScore)
{
#ifdef PANCAKE_MAP_CLR_DEBUG_2
    std::cerr << "(ReChainSeedHits_) Starting to rechain the seed hits.\n";
#endif

    std::vector<std::unique_ptr<ChainedRegion>> newChainedRegions;

    // Merge all remaining seed hits.
    std::vector<SeedHit> hits2;
    for (size_t i = 0; i < chainedRegions.size(); ++i) {
        auto& region = chainedRegions[i];
        if (region->priority > 1) {
            continue;
        }
        hits2.insert(hits2.end(), region->chain.hits.begin(), region->chain.hits.end());
    }

    // Sort the hits by coordinates.
    // IMPORTANT: This needs to sort by target, and if target coords are identical then by query.
    std::sort(hits2.begin(), hits2.end(), [](const SeedHit& a, const SeedHit& b) {
        return std::tuple(a.targetId, a.targetRev, a.targetPos, a.queryPos) <
               std::tuple(b.targetId, b.targetRev, b.targetPos, b.queryPos);
    });

    auto groups = GroupByTargetAndStrand(hits2);

    for (const auto& group : groups) {
#ifdef PANCAKE_MAP_CLR_DEBUG_2
        std::cerr << "(ReChainSeedHits_) group: " << group << ", span = " << group.Span() << "\n";
        for (int32_t i = group.start; i < group.end; ++i) {
            std::cerr << "    [hit i = " << i << "] " << hits2[i] << "\n";
        }
        std::cerr << "\n";
#endif

        // DP Chaining of the filtered hits to remove outliers.
        std::vector<ChainedHits> chains = ChainHits(
            &hits2[group.start], group.end - group.start, chainMaxSkip, chainMaxPredecessors,
            maxGap, chainBandwidth, minNumSeeds, minCoveredBases, minDPScore);

        // Accumulate chains and their mapped regions.
        for (size_t i = 0; i < chains.size(); ++i) {
            const auto& chain = chains[i];
            int32_t numHitsInChain = chain.hits.size();

            // Filter.
            if (numHitsInChain < minNumSeeds) {
                continue;
            }

            // Create a new chained region.
            auto ovl = MakeOverlap_(chain.hits, queryId, queryLen, index, 0, numHitsInChain, 0,
                                    numHitsInChain - 1);
            auto chainedRegion = std::make_unique<ChainedRegion>();
            chainedRegion->chain = std::move(chain);
            chainedRegion->mapping = std::move(ovl);
            newChainedRegions.emplace_back(std::move(chainedRegion));
        }
    }

    return newChainedRegions;
}

std::vector<std::unique_ptr<ChainedRegion>> MapperCLR::ChainAndMakeOverlap_(
    const PacBio::Pancake::SeedIndex& index, const std::vector<SeedHit>& hits,
    const std::vector<PacBio::Pancake::Range>& hitGroups, int32_t queryId, int32_t queryLen,
    int32_t chainMaxSkip, int32_t chainMaxPredecessors, int32_t maxGap, int32_t chainBandwidth,
    int32_t minNumSeeds, int32_t minCoveredBases, int32_t minDPScore, bool useLIS)
{
    // Comparison function to sort the seed hits for LIS.
    // IMPORTANT: This needs to sort by target, and if target coords are identical then by query.
    auto ComparisonSort = [](const SeedHit& a, const SeedHit& b) -> bool {
        return std::pair(a.targetPos, a.queryPos) < std::pair(b.targetPos, b.queryPos);
    };

    // Comparison function to compute the LIS.
    // IMPORTANT: This needs to always return the upper-left element as the smaller one.
    std::function<bool(const SeedHit& a, const SeedHit& b)> ComparisonLIS = [](const SeedHit& a,
                                                                               const SeedHit& b) {
        return (a.queryPos < b.queryPos && a.targetPos < b.targetPos);
    };

    // Process each diagonal bin to get the final chains.
    std::vector<std::unique_ptr<ChainedRegion>> allChainedRegions;
    for (const auto& range : hitGroups) {
#ifdef PANCAKE_MAP_CLR_DEBUG_2
        std::cerr << "(ChainAndMakeOverlap_) [range] start = " << range.start
                  << ", end = " << range.end << ", span = " << range.Span() << "\n";
#endif
        // Skip diagonals with insufficient hits.
        if (range.Span() < minNumSeeds) {
#ifdef PANCAKE_MAP_CLR_DEBUG_2
            std::cerr << "  -> Skipping, range.Span() < minNumSeeds\n";
#endif
            continue;
        }

        // Get a copy of the seed hits so that we can sort without touching the input.
        std::vector<SeedHit> groupHits(hits.begin() + range.start, hits.begin() + range.end);

        // Groups are already groupped by target ID and strand, so only sorting by coordinates is enough.
        // Hits have previously been sorted by diagonals, and not by coordinates, so we need to sort again
        // to get them in proper order.
        std::sort(groupHits.begin(), groupHits.end(), ComparisonSort);

        // Perform chaining.
        std::vector<ChainedHits> chains;

        if (useLIS) {
            // Longest Increasing Subsequence.
            std::vector<PacBio::Pancake::SeedHit> lisHits =
                istl::LIS(groupHits, 0, groupHits.size(), ComparisonLIS);

            // DP Chaining of the filtered hits to remove outliers.
            chains = ChainHits(&lisHits[0], lisHits.size(), chainMaxSkip, chainMaxPredecessors,
                               maxGap, chainBandwidth, minNumSeeds, minCoveredBases, minDPScore);

#ifdef PANCAKE_MAP_CLR_DEBUG_2
            std::cerr << "  - Hits before LIS:\n";
            for (size_t ii = 0; ii < groupHits.size(); ++ii) {
                std::cerr << "    [groupHits ii = " << ii << "] " << groupHits[ii] << "\n";
            }
            std::cerr << "  - using LIS.\n";
            std::cerr << "  - lisHits.size() = " << lisHits.size() << "\n";
            std::cerr << "  - Hits after LIS:\n";
            for (size_t ii = 0; ii < lisHits.size(); ++ii) {
                std::cerr << "    [lisHits ii = " << ii << "] " << lisHits[ii] << "\n";
            }
#endif
        } else {
            chains = ChainHits(&groupHits[0], groupHits.size(), chainMaxSkip, chainMaxPredecessors,
                               maxGap, chainBandwidth, minNumSeeds, minCoveredBases, minDPScore);
#ifdef PANCAKE_MAP_CLR_DEBUG_2
            std::cerr << "  - not using LIS.\n";
#endif
        }

#ifdef PANCAKE_MAP_CLR_DEBUG_2
        std::cerr << "  - Adding chains, chains.size() = " << chains.size() << "\n";
#endif

        // Accumulate chains and their mapped regions.
        for (size_t i = 0; i < chains.size(); ++i) {
#ifdef PANCAKE_MAP_CLR_DEBUG_2
            std::cerr << "    [chain i = " << i << "]\n";
#endif

            const auto& chain = chains[i];
            int32_t numHitsInChain = chain.hits.size();
            // Filter.
            if (numHitsInChain < minNumSeeds) {
#ifdef PANCAKE_MAP_CLR_DEBUG_2
                std::cerr << "      -> Skipping chain, numHitsInChain < minNumSeeds\n";
#endif
                continue;
            }
            // Create a new chained region.
            auto ovl = MakeOverlap_(chain.hits, queryId, queryLen, index, 0, numHitsInChain, 0,
                                    numHitsInChain - 1);
            auto chainedRegion = std::make_unique<ChainedRegion>();
            chainedRegion->chain = std::move(chain);
            chainedRegion->mapping = std::move(ovl);
            allChainedRegions.emplace_back(std::move(chainedRegion));
        }
    }
    return allChainedRegions;
}

OverlapPtr MapperCLR::MakeOverlap_(const std::vector<SeedHit>& sortedHits, int32_t queryId,
                                   int32_t queryLen, const PacBio::Pancake::SeedIndex& index,
                                   int32_t beginId, int32_t endId, int32_t minTargetPosId,
                                   int32_t maxTargetPosId)
{

    const auto& beginHit = sortedHits[minTargetPosId];
    const auto& endHit = sortedHits[maxTargetPosId];

    const int32_t targetId = beginHit.targetId;
    const int32_t numSeeds = endId - beginId;

    if (endHit.targetId != beginHit.targetId) {
        std::ostringstream oss;
        oss << "The targetId of the first and last seed does not match, in MakeOverlap. "
               "beginHit.targetId "
            << beginHit.targetId << ", endHit.targetId = " << endHit.targetId;
        throw std::runtime_error(oss.str());
    }

    const float score = numSeeds;
    const float identity = 0.0;
    const int32_t editDist = -1;

    const int32_t targetLen = index.GetSequenceLength(targetId);

    OverlapPtr ret =
        createOverlap(queryId, targetId, score, identity, beginHit.targetRev, beginHit.queryPos,
                      endHit.queryPos, queryLen, false, beginHit.targetPos, endHit.targetPos,
                      targetLen, editDist, numSeeds, OverlapType::Unknown, OverlapType::Unknown);

    ret->NormalizeStrand();

    return ret;
}

void MapperCLR::LongMergeChains_(std::vector<std::unique_ptr<ChainedRegion>>& chainedRegions,
                                 int32_t maxGap)
{
    if (maxGap < 0) {
        throw std::runtime_error("(LongMergeChains_) maxGap cannot be negative. maxGap = " +
                                 std::to_string(maxGap));
    }

    if (chainedRegions.empty()) {
        return;
    }

    std::vector<int32_t> candidates;
    for (size_t i = 0; i < chainedRegions.size(); ++i) {
        if (chainedRegions[i] == nullptr) {
            continue;
        }
        auto& cr = chainedRegions[i];
        if (cr->priority < 0 || cr->priority > 1) {
            continue;
        }
        candidates.emplace_back(i);
    }

    std::sort(candidates.begin(), candidates.end(), [&](const auto& a, const auto& b) {
        const auto& ar = chainedRegions[a];
        const auto& br = chainedRegions[b];
        return std::tuple(ar->mapping->Bid, ar->mapping->Brev, ar->mapping->Astart,
                          ar->mapping->Bstart) < std::tuple(br->mapping->Bid, br->mapping->Brev,
                                                            br->mapping->Astart,
                                                            br->mapping->Bstart);
    });

    std::unordered_set<int32_t> doSort;

    int32_t lastId = candidates[0];
    for (int32_t i = 1; i < static_cast<int32_t>(candidates.size()); ++i) {
        const int32_t currId = candidates[i];
        auto& last = chainedRegions[lastId];
        const auto& curr = chainedRegions[currId];

        // Skip if there is an overlap (or if the target is wrong), to preserve colinearity.
        if (curr->mapping->Bid != last->mapping->Bid ||
            curr->mapping->Brev != last->mapping->Brev ||
            curr->mapping->Astart <= last->mapping->Aend ||
            curr->mapping->Bstart <= last->mapping->Bend) {
            lastId = currId;
            continue;
        }

        const int32_t gap = std::abs((curr->mapping->Bstart - last->mapping->Bend) -
                                     (curr->mapping->Astart - last->mapping->Aend));

        if (gap > maxGap) {
            lastId = currId;
            continue;
        }

        last->chain.hits.insert(last->chain.hits.end(), curr->chain.hits.begin(),
                                curr->chain.hits.end());
        last->chain.score += curr->chain.score;
        last->chain.coveredBasesQuery += curr->chain.coveredBasesQuery;
        last->chain.coveredBasesTarget += curr->chain.coveredBasesTarget;

        last->mapping->Aend = curr->mapping->Aend;
        last->mapping->Bend = curr->mapping->Bend;
        last->mapping->Score += curr->mapping->Score;
        last->mapping->NumSeeds += curr->mapping->NumSeeds;

        chainedRegions[currId] = nullptr;

        doSort.emplace(lastId);
    }

    // Sort only the extended chains.
    for (const auto& i : doSort) {
        if (chainedRegions[i] == nullptr) {
            continue;
        }
        std::sort(chainedRegions[i]->chain.hits.begin(), chainedRegions[i]->chain.hits.end());
    }

    // Remove the merged ones.
    std::vector<std::unique_ptr<ChainedRegion>> ret;
    for (size_t i = 0; i < chainedRegions.size(); ++i) {
        if (chainedRegions[i] == nullptr) {
            continue;
        }
        if (chainedRegions[i]->mapping == nullptr) {
            continue;
        }
        ret.emplace_back(std::move(chainedRegions[i]));
    }

    std::swap(chainedRegions, ret);
}

}  // namespace Pancake
}  // namespace PacBio
