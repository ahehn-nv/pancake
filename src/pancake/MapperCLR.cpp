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

void DebugPrintChainedRegion(std::ostream& oss, int32_t regionId, const ChainedRegion& cr);

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

std::vector<MapperCLRResult> MapperCLR::Map(const std::vector<std::string>& targetSeqs,
                                            const std::vector<std::string>& querySeqs)
{
    // Collect all seeds for the target sequences.
    std::vector<PacBio::Pancake::Int128t> seeds;
    std::vector<int32_t> sequenceLengths;
    sequenceLengths.reserve(targetSeqs.size());
    for (int32_t recordId = 0; recordId < static_cast<int32_t>(targetSeqs.size()); ++recordId) {
        const auto& record = targetSeqs[recordId];
        const uint8_t* seq = reinterpret_cast<const uint8_t*>(record.data());
        int32_t seqLen = record.size();
        sequenceLengths.emplace_back(seqLen);
        std::vector<PacBio::Pancake::Int128t> newSeeds;
        int rv = SeedDB::GenerateMinimizers(
            newSeeds, seq, seqLen, 0, recordId, settings_.seedParams.KmerSize,
            settings_.seedParams.MinimizerWindow, settings_.seedParams.Spacing,
            settings_.seedParams.UseRC, settings_.seedParams.UseHPCForSeedsOnly,
            settings_.seedParams.MaxHPCLen);
        if (rv)
            throw std::runtime_error("Generating minimizers failed for the target sequence, id = " +
                                     std::to_string(recordId));
        seeds.insert(seeds.end(), newSeeds.begin(), newSeeds.end());
    }

    // Construct the index.
    SeedIndex seedIndex(settings_.seedParams, sequenceLengths, std::move(seeds));

    // Calculate the seed frequency statistics, needed for the cutoff.
    int64_t freqMax = 0;
    int64_t freqCutoff = 0;
    double freqAvg = 0.0;
    double freqMedian = 0.0;
    seedIndex.ComputeFrequencyStats(settings_.freqPercentile, freqMax, freqAvg, freqMedian,
                                    freqCutoff);

    // Run mapping for each query.
    std::vector<MapperCLRResult> results;
    for (int32_t queryId = 0; queryId < static_cast<int32_t>(querySeqs.size()); ++queryId) {
        const auto& query = querySeqs[queryId];

#ifdef PANCAKE_MAP_CLR_DEBUG
        PBLOG_INFO << "New query: queryId = " << queryId << ", length = " << query.size();
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
            throw std::runtime_error("Generating minimizers failed for the query sequence, id = " +
                                     std::to_string(queryId));

        auto queryResults = Map(targetSeqs, seedIndex, query, querySeeds, queryId, freqCutoff);
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
    oss << "[regionId " << regionId << "] size = " << cr.chain.hits.size()
        << ", score = " << cr.chain.score << ", covQ = " << cr.chain.coveredBasesQuery
        << ", covT = " << cr.chain.coveredBasesTarget << ", priority = " << cr.priority
        << ", isSuppl = " << (cr.isSupplementary ? "true" : "false");
    // int32_t numHitsInChain = chain.hits.size();
    // auto ovl = MakeOverlap_(chain.hits, 0, querySeq.size(), index, 0, numHitsInChain, 0,
    //                         numHitsInChain - 1);
    oss << " " << OverlapWriterBase::PrintOverlapAsM4(cr.mapping, "", "", true, false);
    oss << ", diagStart = " << (cr.mapping->Astart - cr.mapping->Bstart);
    oss << ", diagEnd = " << (cr.mapping->Aend - cr.mapping->Bend);
}

MapperCLRResult MapperCLR::Map(const std::vector<std::string>& targetSeqs,
                               const PacBio::Pancake::SeedIndex& index, const std::string& querySeq,
                               const std::vector<PacBio::Pancake::Int128t>& querySeeds,
                               const int32_t queryId, int64_t freqCutoff)
{
    // #ifdef PANCAKE_DEBUG
    //     PBLOG_INFO << "Mapping query ID = " << querySeq.Id() << ", header = " << querySeq.Name();
    // #endif
    if (static_cast<int64_t>(querySeq.size()) < settings_.minQueryLen) {
        return {};
    }

    // Reverse the query sequence and make a helper array.
    const std::string querySeqRev =
        PacBio::Pancake::ReverseComplement(querySeq, 0, querySeq.size());
    const std::array<const char*, 2> queryArray = {&querySeq[0], &querySeqRev[0]};
    const int32_t queryLen = querySeq.size();

    TicToc ttCollectHits;
    std::vector<SeedHit> hits;
    index.CollectHits(&querySeeds[0], querySeeds.size(), queryLen, hits, freqCutoff);
    ttCollectHits.Stop();

    TicToc ttSortHits;
    std::sort(hits.begin(), hits.end(), [](const auto& a, const auto& b) {
        return PackSeedHitWithDiagonalToTuple(a) < PackSeedHitWithDiagonalToTuple(b);
    });
    ttSortHits.Stop();

    // PBLOG_INFO << "Hits: " << hits.size();

    TicToc ttDiagonalGroup;
    auto groups = DiagonalGroup_(hits, settings_.chainBandwidth, true);
    ttDiagonalGroup.Stop();

#ifdef PANCAKE_MAP_CLR_DEBUG_2
    std::cerr << "Diagonal groups:\n";
    for (size_t i = 0; i < groups.size(); ++i) {
        const int32_t firstDiag = hits[groups[i].start].Diagonal();
        const int32_t lastDiag = hits[groups[i].end - 1].Diagonal();
        std::cerr << "[queryId = " << queryId << ", group " << i << "] start = " << groups[i].start
                  << ", end = " << groups[i].end << ", diagStart = " << firstDiag
                  << ", diagEnd = " << lastDiag << "\n";
    }
#endif

#ifdef PANCAKE_WRITE_SCATTERPLOT
    for (size_t i = 0; i < groups.size(); ++i) {
        const auto& g = groups[i];
        const int32_t targetId = hits[g.start].targetId;
        const int32_t targetLen = targetSeqs[targetId].size();
        WriteSeedHits_("temp-debug/hits-q" + std::to_string(queryId) + "-0-diag.csv", hits, g.start,
                       g.end, i, "query" + std::to_string(queryId), queryLen,
                       "target" + std::to_string(targetId), targetLen, (i > 0));
    }
#endif

    // Process each diagonal bin to get the final chains.
    TicToc ttChain;
    std::vector<std::unique_ptr<ChainedRegion>> allChainedRegions = ChainAndMakeOverlap_(
        index, hits, groups, queryId, queryLen, settings_.chainMaxSkip,
        settings_.chainMaxPredecessors, settings_.maxGap, settings_.chainBandwidth,
        settings_.minNumSeeds, settings_.minCoveredBases, settings_.minDPScore, settings_.useLIS);
    ttChain.Stop();

#ifdef PANCAKE_MAP_CLR_DEBUG_2
    std::cerr << "After ChainAndMakeOverlap_: allChainedRegions.size() = "
              << allChainedRegions.size() << "\n";
#endif

#ifdef PANCAKE_WRITE_SCATTERPLOT
    for (size_t i = 0; i < allChainedRegions.size(); ++i) {
        auto& region = allChainedRegions[i];
        WriteSeedHits_(
            "temp-debug/hits-q" + std::to_string(queryId) + "-1-chained.csv", region->chain.hits, 0,
            region->chain.hits.size(), i, "query" + std::to_string(queryId), region->mapping->Alen,
            "target" + std::to_string(region->mapping->Bid), region->mapping->Blen, (i > 0));
    }
#endif

    // Take the remaining regions, merge all seed hits, and rechain.
    // Needed because diagonal chaining was greedy and a wide window could have split
    // otherwise good chains. On the other hand, diagonal binning was needed for speed
    // in low-complexity regions.
    allChainedRegions =
        ReChainSeedHits_(allChainedRegions, index, queryId, queryLen, settings_.chainMaxSkip,
                         settings_.chainMaxPredecessors, settings_.maxGap, settings_.chainBandwidth,
                         settings_.minNumSeeds, settings_.minCoveredBases, settings_.minDPScore);

#ifdef PANCAKE_WRITE_SCATTERPLOT
    for (size_t i = 0; i < allChainedRegions.size(); ++i) {
        auto& region = allChainedRegions[i];
        if (region->priority > 1) {
            continue;
            // break;
        }
        WriteSeedHits_("temp-debug/hits-q" + std::to_string(queryId) + "-2-rechained.csv",
                       region->chain.hits, 0, region->chain.hits.size(), i,
                       "query" + std::to_string(queryId), region->mapping->Alen,
                       "target" + std::to_string(region->mapping->Bid), region->mapping->Blen,
                       (i > 0));
    }
#endif

    // Sort all chains in descending order of the number of hits.
    std::sort(allChainedRegions.begin(), allChainedRegions.end(),
              [](const auto& a, const auto& b) { return a->chain.score > b->chain.score; });

    // Secondary/supplementary flagging.
    WrapFlagSecondaryAndSupplementary_(allChainedRegions, settings_.secondaryAllowedOverlapFraction,
                                       settings_.secondaryMinScoreFraction);

#ifdef PANCAKE_WRITE_SCATTERPLOT
    for (size_t i = 0; i < allChainedRegions.size(); ++i) {
        auto& region = allChainedRegions[i];
        if (region->priority > 1) {
            continue;
            // break;
        }
        WriteSeedHits_("temp-debug/hits-q" + std::to_string(queryId) + "-3-secondary.csv",
                       region->chain.hits, 0, region->chain.hits.size(), i,
                       "query" + std::to_string(queryId), region->mapping->Alen,
                       "target" + std::to_string(region->mapping->Bid), region->mapping->Blen,
                       (i > 0));
    }
#endif

    // Merge long gaps.
    LongMergeChains_(allChainedRegions, settings_.maxGap);

    // Again relabel, because some chains are longer now.
    WrapFlagSecondaryAndSupplementary_(allChainedRegions, settings_.secondaryAllowedOverlapFraction,
                                       settings_.secondaryMinScoreFraction);

    // Sort all chains by priority and then score.
    std::sort(allChainedRegions.begin(), allChainedRegions.end(),
              [](const std::unique_ptr<ChainedRegion>& a, const std::unique_ptr<ChainedRegion>& b) {
                  auto at = std::tuple<int32_t, bool, int32_t>(a->priority, a->isSupplementary,
                                                               a->chain.score);
                  auto bt = std::tuple<int32_t, bool, int32_t>(b->priority, b->isSupplementary,
                                                               b->chain.score);
                  return at < bt;
              });

    // Filter seed hits.
    for (size_t i = 0; i < allChainedRegions.size(); ++i) {
        auto& region = allChainedRegions[i];
        if (region->priority > 1) {
            continue;
            // break;
        }
        ChainedHits newChain = RefineBadEnds(
            region->chain, settings_.alnParamsGlobal.alignBandwidth, settings_.minDPScore * 2);
        // auto newChain = region->chain;
        newChain = RefineChainedHits(newChain, 10, 40, settings_.maxGap / 2, 10);
        newChain = RefineChainedHits2(newChain, 30, settings_.maxGap / 2);

#ifdef PANCAKE_MAP_CLR_DEBUG_2
        const double divergenceOld = ComputeChainDivergence(region->chain.hits);
        const double divergenceNew = ComputeChainDivergence(newChain.hits);
        const int32_t numHitsOld = region->chain.hits.size();
        const int32_t numHitsNew = newChain.hits.size();
#endif

        std::swap(region->chain, newChain);

#ifdef PANCAKE_MAP_CLR_DEBUG_2
        std::cerr << "(refined i = " << i << " / " << allChainedRegions.size() << ") ";
        DebugPrintChainedRegion(std::cerr, i, *region);
        std::cerr << ", divergOld = " << divergenceOld << ", divergNew = " << divergenceNew
                  << ", hitsOld = " << numHitsOld << ", hitsNew = " << numHitsNew << "\n";
#endif
    }

#ifdef PANCAKE_WRITE_SCATTERPLOT
    for (size_t i = 0; i < allChainedRegions.size(); ++i) {
        auto& region = allChainedRegions[i];
        if (region->priority > 1) {
            continue;
            // break;
        }
        WriteSeedHits_(
            "temp-debug/hits-q" + std::to_string(queryId) + "-4-refined.csv", region->chain.hits, 0,
            region->chain.hits.size(), i, "query" + std::to_string(queryId), region->mapping->Alen,
            "target" + std::to_string(region->mapping->Bid), region->mapping->Blen, (i > 0));
    }
#endif

#ifdef PANCAKE_MAP_CLR_DEBUG

    std::cerr << "hits.size() = " << hits.size() << "\n";

    // std::sort(hits.begin(), hits.end(),
    //           [](const auto& a, const auto& b) { return a.PackTo128() < b.PackTo128(); });
    // for (size_t j = 0; j < hits.size(); ++j) {
    //     const auto& hit = hits[j];
    //     std::cerr << "[hit " << j << "] " << hit << "\n";
    // }

    // std::cerr << "groups.size() = " << groups.size() << "\n";
    // for (size_t i = 0; i < groups.size(); ++i) {
    //     std::cerr << "   [group " << i << "] size = " << groups[i].Span() << "\n";
    // }

    std::cerr << "allChainedRegions before skipping secondary, allChainedRegions.size() = "
              << allChainedRegions.size() << ":\n";
    for (size_t i = 0; i < allChainedRegions.size(); ++i) {
        std::cerr << "(final all) ";
        DebugPrintChainedRegion(std::cerr, i, *allChainedRegions[i]);
        std::cerr << "\n";
    }

// for (size_t i = 0; i < allChainedRegions.size(); ++i) {
//     std::cerr << "[queryId = " << queryId << ", hits for chainedRegion " << i << "] score = " << allChainedRegions[i]->chain.score << ", covBasesQuery = " << allChainedRegions[i]->chain.coveredBasesQuery << ", covBasesTarget = " << allChainedRegions[i]->chain.coveredBasesTarget << "\n";
//     for (size_t j = 0; j < allChainedRegions[i]->chain.hits.size(); ++j) {
//         const auto& hit = allChainedRegions[i]->chain.hits[j];
//         std::cerr << "    [hit " << j << "] tid = " << hit.targetId << ", trev = " << hit.targetRev << ", tpos = " << hit.targetPos << ", qpos = " << hit.queryPos << "\n";
//     }
// }
#endif

    // Filter out the mappings.
    MapperCLRResult result;
    int32_t numSelectedSecondary = 0;
    for (size_t i = 0; i < allChainedRegions.size(); ++i) {
        const auto& region = allChainedRegions[i];
        if (region->priority > 1) {
            continue;
            // break;
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

#ifdef PANCAKE_MAP_CLR_DEBUG
    std::cerr << "Results:\n";
    for (size_t i = 0; i < result.mappings.size(); ++i) {
        std::cerr << "(results) ";
        DebugPrintChainedRegion(std::cerr, i, *result.mappings[i]);
        std::cerr << "\n";
    }

    std::cerr << "Results with hits:\n";
    for (size_t i = 0; i < result.mappings.size(); ++i) {
        std::cerr << "(result hits) ";
        DebugPrintChainedRegion(std::cerr, i, *result.mappings[i]);
        std::cerr << "\n";
        const auto& m = result.mappings[i];
        std::cerr << "[queryId = " << queryId << ", hits for result " << i
                  << "] score = " << m->chain.score
                  << ", covBasesQuery = " << m->chain.coveredBasesQuery
                  << ", covBasesTarget = " << m->chain.coveredBasesTarget << "\n";
        for (size_t j = 0; j < m->chain.hits.size(); ++j) {
            const auto& hit = m->chain.hits[j];
            std::cerr << "    [hit " << j << "] tid = " << hit.targetId
                      << ", trev = " << hit.targetRev << ", tpos = " << hit.targetPos
                      << ", qpos = " << hit.queryPos << "\n";
        }
    }
#endif

    if (settings_.align) {
#ifdef PANCAKE_MAP_CLR_DEBUG_ALIGN
        std::cerr << "Aligning.\n";
        std::cerr << "alignerTypeGlobal = " << AlignerTypeToString(settings_.alignerTypeGlobal)
                  << "\n";
        std::cerr << "alignerTypeExt = " << AlignerTypeToString(settings_.alignerTypeExt) << "\n";
#endif

        TicToc ttAlign;

        for (size_t i = 0; i < result.mappings.size(); ++i) {
            // const auto& chain = result.mappings[i]->chain;
            auto& chain = result.mappings[i]->chain;
            auto& ovl = result.mappings[i]->mapping;
            const auto& tSeqFwd = targetSeqs[ovl->Bid];

#ifdef PANCAKE_MAP_CLR_DEBUG_ALIGN
            std::cerr << "[Aligning overlap] "
                      << OverlapWriterBase::PrintOverlapAsM4(ovl, "", "", true, false) << "\n";
#endif

            // auto newHits =
            //     FilterSeedHitsByExtension_(tSeq, querySeq, chain.hits, 50, 2, 30, sesScratch_);
            // std::swap(chain.hits, newHits);
            // chain = RefineChainedHits(chain, k, 10, 40, settings_.maxGap / 2, 10);
            // chain = RefineChainedHits2(chain, k, 30, settings_.maxGap / 2);

            // Use a custom aligner to align.
            TicToc ttAlignBetweenSeeds;
            ovl =
                AlignOverlapSeeded(ovl, chain.hits, tSeqFwd.c_str(), tSeqFwd.size(), queryArray[0],
                                   queryArray[1], queryLen, settings_.minAlignmentSpan,
                                   settings_.maxFlankExtensionDist, alignerGlobal_, alignerExt_);
            ttAlignBetweenSeeds.Stop();

#ifdef PANCAKE_MAP_CLR_DEBUG_ALIGN
            std::cerr << "ttAlignBetweenSeeds = " << ttAlignBetweenSeeds.GetCpuMillisecs()
                      << " ms\n";

            if (ovl != nullptr) {
                Alignment::DiffCounts diffs = CigarDiffCounts(ovl->Cigar);
                float idt = 0.0;
                int32_t editDist = 0;
                diffs.Identity(false, false, idt, editDist);
                std::cerr << "Diffs: " << diffs << ", editDist: " << editDist
                          << ", qspan = " << (diffs.numEq + diffs.numX + diffs.numI)
                          << ", tspan = " << (diffs.numEq + diffs.numX + diffs.numD)
                          << ", identity = " << idt * 100.0f << "% \n";
                std::cerr << "CIGAR: " << ovl->Cigar.ToStdString() << "\n";
            } else {
                std::cerr << "ovl == nullptr\n";
            }
#endif
            // TicToc ttAlignGlobal;
            // std::vector<PacBio::Pancake::SeedHit> dummyHits;
            // dummyHits.emplace_back(chain.hits.front());
            // dummyHits.emplace_back(chain.hits.back());
            // AlignBetweenSeeds_(tSeq, querySeq, dummyHits, settings_.alignmentMaxD,
            //                 settings_.alignmentBandwidth, sesScratch_);
            // ttAlignGlobal.Stop();
            // std::cerr << "ttAlignGlobal = " << ttAlignGlobal.GetCpuMillisecs() << " ms\n";
            // break;
        }

        ttAlign.Stop();
#ifdef PANCAKE_MAP_CLR_DEBUG_ALIGN
        std::cerr << "ttAlign = " << ttAlign.GetCpuMillisecs() << " ms\n";
#endif
    }

#ifdef PANCAKE_MAP_CLR_DEBUG_ALIGN
    std::cerr << "\n";
#endif

    // Remove any mappings which didn't align.
    size_t numValid = 0;
    for (size_t i = 0; i < result.mappings.size(); ++i) {
        if (result.mappings[i]->mapping == nullptr) {
            continue;
        }
        if (i != numValid) {
            std::swap(result.mappings[i], result.mappings[numValid]);
        }
        ++numValid;
    }
    result.mappings.resize(numValid);

#ifdef PANCAKE_WRITE_SCATTERPLOT
    for (size_t i = 0; i < result.mappings.size(); ++i) {
        auto& region = result.mappings[i];
        if (region->mapping == nullptr) {
            continue;
        }
        if (region->priority > 1) {
            continue;
        }
        WriteSeedHits_(
            "temp-debug/hits-q" + std::to_string(queryId) + "-5-results.csv", region->chain.hits, 0,
            region->chain.hits.size(), i, "query" + std::to_string(queryId), region->mapping->Alen,
            "target" + std::to_string(region->mapping->Bid), region->mapping->Blen, (i > 0));
    }
#endif

    return result;
}

void MapperCLR::WrapFlagSecondaryAndSupplementary_(
    std::vector<std::unique_ptr<ChainedRegion>>& allChainedRegions,
    double secondaryAllowedOverlapFraction, double secondaryMinScoreFraction)
{
    /*
     * Edits the objects in place.
    */
    // Copy the overlaps so we can satisfy the FlagSecondaryAndSupplementary API.
    std::vector<OverlapPtr> tmpOverlaps;
    for (size_t i = 0; i < allChainedRegions.size(); ++i) {
        tmpOverlaps.emplace_back(createOverlap(allChainedRegions[i]->mapping));
    }
    // Flag the secondary and supplementary overlaps.
    std::vector<OverlapPriority> overlapPriorities = FlagSecondaryAndSupplementary(
        tmpOverlaps, secondaryAllowedOverlapFraction, secondaryMinScoreFraction);
    // Set the flags.
    for (size_t i = 0; i < allChainedRegions.size(); ++i) {
        auto& cr = allChainedRegions[i];
        cr->mapping->IsSecondary = (overlapPriorities[i].priority > 0);
        cr->mapping->IsSupplementary = overlapPriorities[i].isSupplementary;
        cr->priority = overlapPriorities[i].priority;
        cr->isSupplementary = overlapPriorities[i].isSupplementary;
    }
}

void MapperCLR::WriteSeedHits_(const std::string& outPath, const std::vector<SeedHit>& hits,
                               size_t hitsStart, size_t hitsEnd, int32_t hitsId,
                               const std::string& queryName, int64_t queryLength,
                               const std::string& targetName, int64_t targetLength, bool append)
{
    if (hitsStart >= hits.size() || hitsEnd > hits.size() || hitsStart > hitsEnd) {
        std::ostringstream oss;
        oss << "Invalid hitsStart and/or hitsEnd! hitsStart = " << hitsStart
            << ", hitsEnd = " << hitsEnd << ", hits.size() = " << hits.size();
        throw std::runtime_error(oss.str());
    }

    std::ofstream ofs;
    if (append) {
        ofs = std::ofstream(outPath, std::ios::app);
    } else {
        ofs = std::ofstream(outPath);
    }
    if (ofs.is_open() == false) {
        // Don't throw. This is a hidden feature which will only work if a user knows which folder to create.
        return;
        // throw std::runtime_error("Could not open file '" + outPath +
        //                          "' for writing! In MapperCLR::WriteSeedHits.");
    }
    if (append == false) {
        ofs << queryName.c_str() << "\t" << queryLength << "\t" << targetName.c_str() << "\t"
            << targetLength << "\n";
    }
    for (size_t j = hitsStart; j < hitsEnd; ++j) {
        ofs << hits[j].queryPos << "\t" << hits[j].targetPos << "\t" << hitsId << "\n";
        // ofs << (hits[j].queryPos + kmerSize) << "\t" << (hits[j].targetPos + kmerSize) << "\t"
        //     << hitsId << "\n";
    }
}

std::vector<std::unique_ptr<ChainedRegion>> MapperCLR::ReChainSeedHits_(
    const std::vector<std::unique_ptr<ChainedRegion>>& chainedRegions,
    const PacBio::Pancake::SeedIndex& index, int32_t queryId, int32_t queryLen,
    int32_t chainMaxSkip, int32_t chainMaxPredecessors, int32_t maxGap, int32_t chainBandwidth,
    int32_t minNumSeeds, int32_t minCoveredBases, int32_t minDPScore)
{
    std::vector<std::unique_ptr<ChainedRegion>> newChainedRegions;

    // Collect all remaining seed hits.
    std::vector<SeedHit> hits2;
    for (size_t i = 0; i < chainedRegions.size(); ++i) {
        auto& region = chainedRegions[i];
        if (region->priority > 1) {
            continue;
            // break;
        }
        hits2.insert(hits2.end(), region->chain.hits.begin(), region->chain.hits.end());
    }

    // Sort the hits by coordinates.
    std::sort(hits2.begin(), hits2.end(), [](const SeedHit& a, const SeedHit& b) {
        // return a.queryPos < b.queryPos || (a.queryPos == b.queryPos && a.targetPos < b.targetPos);
        return std::tuple(a.targetId, a.targetRev, a.targetPos, a.queryPos) <
               std::tuple(b.targetId, b.targetRev, b.targetPos, b.queryPos);
        // return std::tuple(a.targetId, a.targetRev, a.queryPos, a.targetPos) <
        //        std::tuple(b.targetId, b.targetRev, b.queryPos, b.targetPos);

        // return a.targetPos < b.targetPos || (a.targetPos == b.targetPos && a.queryPos < b.queryPos);
    });

    auto groups = GroupByTargetAndStrand_(hits2);

#ifdef PANCAKE_WRITE_SCATTERPLOT
    // {
    //     const int32_t targetId = hits2[0].targetId;
    //     const int32_t targetLen = index.GetSequenceLength(targetId);
    //     WriteSeedHits_("temp-debug/hits-q" + std::to_string(queryId) + "-1.1-hits2.csv", hits2, 0,
    //                 hits2.size(), 0, "query" + std::to_string(queryId), queryLen,
    //                 "target" + std::to_string(targetId), targetLen, 0);
    // }
    for (size_t i = 0; i < groups.size(); ++i) {
        const auto& g = groups[i];
        const int32_t targetId = hits2[g.start].targetId;
        const int32_t targetLen = index.GetSequenceLength(targetId);
        WriteSeedHits_("temp-debug/hits-q" + std::to_string(queryId) + "-1.2-hits2-groups.csv",
                       hits2, g.start, g.end, i, "query" + std::to_string(queryId), queryLen,
                       "target" + std::to_string(targetId), targetLen, (i > 0));
    }
#endif

    for (const auto& group : groups) {
        // DP Chaining of the filtered hits to remove outliers.
        // std::cerr << "group.start = " << group.start << ", group.end = " << group.end << "\n";

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
    auto ComparisonSort = [](const SeedHit& a, const SeedHit& b) -> bool {
        // auto aQpos = a.queryPos;
        // auto aTpos = a.targetPos;
        // auto bQpos = b.queryPos;
        // auto bTpos = b.targetPos;
        // return (aQpos < bQpos || (aQpos == bQpos && aTpos < bTpos));
        return std::pair(a.targetPos, a.queryPos) < std::pair(b.targetPos, b.queryPos);
        // return std::pair(a.queryPos, a.targetPos) < std::pair(b.queryPos, b.targetPos);
    };
    std::function<bool(const SeedHit& a, const SeedHit& b)> ComparisonLIS = [](const SeedHit& a,
                                                                               const SeedHit& b) {
        // This needs to always return the upper-left element as the smaller one.
        return (a.queryPos < b.queryPos && a.targetPos < b.targetPos);
    };

    // Process each diagonal bins to get the final chains.
    std::vector<std::unique_ptr<ChainedRegion>> allChainedRegions;
    for (const auto& range : hitGroups) {
#ifdef PANCAKE_MAP_CLR_DEBUG_2
        std::cerr << "[range] start = " << range.start << ", end = " << range.end
                  << ", span = " << range.Span() << "\n";
#endif
        // Skip diagonals with insufficient hits.
        if (range.Span() < minNumSeeds) {
#ifdef PANCAKE_MAP_CLR_DEBUG_2
            std::cerr << "  -> Skipping, range.Span() < minNumSeeds\n";
#endif
            continue;
        }

        // Longest Increasing Subsequence of the diagonal bin.
        std::vector<SeedHit> groupHits(hits.begin() + range.start, hits.begin() + range.end);
        // Groups are already groupped by target ID and strand, so only sorting by coordinates is enough.
        std::sort(groupHits.begin(), groupHits.end(), ComparisonSort);

        // Perform chaining.
        std::vector<ChainedHits> chains;
        if (useLIS) {
            // Longest Increasing Subsequence.
            std::vector<PacBio::Pancake::SeedHit> lisHits =
                istl::LIS(groupHits, 0, groupHits.size(), ComparisonLIS);

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

            // DP Chaining of the filtered hits to remove outliers.
            chains = ChainHits(&lisHits[0], lisHits.size(), chainMaxSkip, chainMaxPredecessors,
                               maxGap, chainBandwidth, minNumSeeds, minCoveredBases, minDPScore);

#ifdef PANCAKE_MAP_CLR_DEBUG_2
            std::cerr << "  - chains.size() = " << chains.size() << "\n";
#endif
        } else {
            chains = ChainHits(&groupHits[0], groupHits.size(), chainMaxSkip, chainMaxPredecessors,
                               maxGap, chainBandwidth, minNumSeeds, minCoveredBases, minDPScore);
#ifdef PANCAKE_MAP_CLR_DEBUG_2
            std::cerr << "  - not using LIS.\n";
            std::cerr << "  - chains.size() = " << chains.size() << "\n";
#endif
        }

#ifdef PANCAKE_MAP_CLR_DEBUG_2
        std::cerr << "  - Adding chains.\n";
#endif

        // Accumulate chains and their mapped regions.
        for (size_t i = 0; i < chains.size(); ++i) {
#ifdef PANCAKE_MAP_CLR_DEBUG_2
            std::cerr << "    [i = " << i << "]\n";
#endif

            const auto& chain = chains[i];
            int32_t numHitsInChain = chain.hits.size();
            // Filter.
            if (numHitsInChain < minNumSeeds) {
#ifdef PANCAKE_MAP_CLR_DEBUG_2
                std::cerr << "      -> Skipping, numHitsInChain < minNumSeeds\n";
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

    // return createOverlap(queryId, targetId, score, identity, 0, beginHit.queryPos, endHit.queryPos,
    //                      queryLen, beginHit.targetRev, beginHit.targetPos, endHit.targetPos,
    //                      targetLen, editDist, numSeeds, OverlapType::Unknown, OverlapType::Unknown);

    OverlapPtr ret =
        createOverlap(queryId, targetId, score, identity, beginHit.targetRev, beginHit.queryPos,
                      endHit.queryPos, queryLen, false, beginHit.targetPos, endHit.targetPos,
                      targetLen, editDist, numSeeds, OverlapType::Unknown, OverlapType::Unknown);

    ret->NormalizeStrand();

    return ret;
}

std::vector<Range> MapperCLR::GroupByTargetAndStrand_(const std::vector<SeedHit>& sortedHits)
{
    if (sortedHits.empty()) {
        return {};
    }
    const int32_t numHits = static_cast<int32_t>(sortedHits.size());
    int32_t beginId = 0;
    std::vector<Range> groups;
    for (int32_t i = 0; i < numHits; ++i) {
        const auto& prevHit = sortedHits[beginId];
        const auto& currHit = sortedHits[i];
        if (currHit.targetId != prevHit.targetId || currHit.targetRev != prevHit.targetRev) {
            groups.emplace_back(Range{beginId, i});
            beginId = i;
        }
    }
    if ((numHits - beginId) > 0) {
        groups.emplace_back(Range{beginId, numHits});
    }
    return groups;
}

std::vector<Range> MapperCLR::DiagonalGroup_(const std::vector<SeedHit>& sortedHits,
                                             int32_t chainBandwidth, bool overlappingWindows)
{
    if (sortedHits.empty()) {
        return {};
    }

    std::vector<Range> groups;

    const int32_t numHits = static_cast<int32_t>(sortedHits.size());
    int32_t beginId = 0;
    int32_t beginDiag = sortedHits[beginId].Diagonal();

    // This is a combination of <targetPos, queryPos>, intended for simple comparison
    // without defining a custom complicated comparison operator.
    uint64_t minTargetQueryPosCombo = (static_cast<uint64_t>(sortedHits[beginId].targetPos) << 32) |
                                      (static_cast<uint64_t>(sortedHits[beginId].queryPos));
    uint64_t maxTargetQueryPosCombo = minTargetQueryPosCombo;
    // int23_t halfBandwidth = chainBandwidth / 2;

    int32_t firstInBandwidth = 0;

    for (int32_t i = 0; i < numHits; ++i) {
        const auto& beginHit = sortedHits[beginId];
        const auto& currHit = sortedHits[i];
        const int32_t currDiag = currHit.Diagonal();
        const int32_t diagDiff = abs(currDiag - beginDiag);
        const uint64_t targetQueryPosCombo =
            (static_cast<uint64_t>(sortedHits[i].targetPos) << 32) |
            (static_cast<uint64_t>(sortedHits[i].queryPos));

#ifdef PANCAKE_MAP_CLR_DEBUG
        std::cerr << "(dg) [hit " << i << "] tid = " << currHit.targetId
                  << ", trev = " << currHit.targetRev << ", tpos = " << currHit.targetPos
                  << ", qpos = " << currHit.queryPos << ", flag = " << currHit.flags
                  << ", diag = " << currDiag << ", firstInBandwidth = " << firstInBandwidth << "\n";
#endif

        if (currHit.targetId != beginHit.targetId || currHit.targetRev != beginHit.targetRev ||
            diagDiff > chainBandwidth) {

            if (overlappingWindows) {
                groups.emplace_back(Range{firstInBandwidth, i});
            } else {
                groups.emplace_back(Range{beginId, i});
            }

            beginId = i;
            beginDiag = currDiag;

            minTargetQueryPosCombo = maxTargetQueryPosCombo = targetQueryPosCombo;

            // Find the earliest hit which is within the bandwidth window from the current hit.
            if (overlappingWindows) {
                for (; firstInBandwidth < i; ++firstInBandwidth) {
                    const auto& firstHit = sortedHits[firstInBandwidth];
                    const int32_t firstDiag = firstHit.Diagonal();
                    const int32_t diagDiffToFirst = abs(firstDiag - beginDiag);
                    if (currHit.targetId != firstHit.targetId ||
                        currHit.targetRev != firstHit.targetRev ||
                        diagDiffToFirst > chainBandwidth) {
                        continue;
                    }
                    break;
                }
            }

#ifdef PANCAKE_MAP_CLR_DEBUG
            std::cerr << "---\n";
#endif
        }

        // Track the minimum and maximum target positions for each diagonal.
        if (targetQueryPosCombo < minTargetQueryPosCombo) {
            minTargetQueryPosCombo = targetQueryPosCombo;
        }
        if (targetQueryPosCombo > maxTargetQueryPosCombo) {
            maxTargetQueryPosCombo = targetQueryPosCombo;
        }
    }

    if ((numHits - beginId) > 0) {
        groups.emplace_back(Range{beginId, numHits});
    }

    // #ifdef PANCAKE_MAP_CLR_DEBUG
    //     for (const auto& range : groups) {
    //         auto ovl = MakeOverlap_(sortedHits, 0, queryLen, index, range.start, range.end, minPosId,
    //                                 maxPosId);
    //         std::cerr << "(dg) ovl->NumSeeds = " << ovl->NumSeeds
    //                   << ", startI = " << beginId << ", endI = " << numHits
    //                   << ", ovl->ASpan() = " << ovl->ASpan()
    //                   << ", ovl->BSpan() = " << ovl->BSpan()
    //                   << ", ovl->Aid = " << ovl->Aid << ", ovl->Bid = " << ovl->Bid
    //                   << ", beginDiag = " << beginDiag << ", endDiag = " << sortedHits.back().Diagonal()
    //                   << "\n";
    //         std::cerr << "(dg) " << OverlapWriterBase::PrintOverlapAsM4(ovl, "", "", true, false)
    //                   << "\n";
    //     }
    // #endif

    return groups;
}

// std::vector<std::unique_ptr<ChainedRegion>>
void MapperCLR::LongMergeChains_(std::vector<std::unique_ptr<ChainedRegion>>& chainedRegions,
                                 int32_t maxGap)
{
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

    // std::cerr << "candidates.size() = " << candidates.size() << "\n";

    // std::vector<std::unique_ptr<ChainedRegion>> filtered;

    std::unordered_set<int32_t> doSort;

    int32_t lastId = candidates[0];
    for (int32_t i = 1; i < static_cast<int32_t>(candidates.size()); ++i) {
        const int32_t currId = candidates[i];
        auto& last = chainedRegions[lastId];
        const auto& curr = chainedRegions[currId];

        // std::cerr << "[i = " << i << "] lastId = " << lastId << ", currId = " << currId << "\n";
        // std::cerr << "    last: " << OverlapWriterBase::PrintOverlapAsM4(last->mapping, "", "", true, false) << "\n";
        // std::cerr << "    curr: " << OverlapWriterBase::PrintOverlapAsM4(curr->mapping, "", "", true, false) << "\n";

        // Skip if there is an overlap (or if the target is wrong), to preserve colinearity.
        if (curr->mapping->Bid != last->mapping->Bid ||
            curr->mapping->Brev != last->mapping->Brev ||
            curr->mapping->Astart < last->mapping->Aend ||
            curr->mapping->Bstart < last->mapping->Bend) {
            lastId = currId;
            // std::cerr << "    -> continue on colinearity\n";
            continue;
        }

        const int32_t gap = std::abs((curr->mapping->Bstart - last->mapping->Bend) -
                                     (curr->mapping->Astart - last->mapping->Aend));

        if (gap > maxGap) {
            lastId = currId;
            // std::cerr << "    -> continue on maxGap\n";
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

        // std::cerr << "    -> merged\n";
    }

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
