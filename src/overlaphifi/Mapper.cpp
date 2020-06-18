// Authors: Ivan Sovic

#include <lib/kxsort/kxsort.h>
#include <pacbio/alignment/AlignmentTools.h>
#include <pacbio/alignment/DiffCounts.h>
#include <pacbio/alignment/SesDistanceBanded.h>
#include <pacbio/overlaphifi/Mapper.h>
#include <pacbio/overlaphifi/OverlapWriterBase.h>
#include <pacbio/seqdb/Util.h>
#include <pacbio/util/RunLengthEncoding.h>
#include <pacbio/util/TicToc.h>
#include <pbcopper/logging/Logging.h>
#include <pbcopper/third-party/edlib.h>
#include <algorithm>
#include <iostream>
#include <pacbio/alignment/Ses2AlignBanded.hpp>
#include <pacbio/alignment/Ses2DistanceBanded.hpp>
#include <pacbio/alignment/SesAlignBanded.hpp>
#include <sstream>

namespace PacBio {
namespace Pancake {

// #define PANCAKE_DEBUG
// #define PANCAKE_DEBUG_ALN

auto AlignWithTraceback(const char* query, size_t queryLen, const char* target, size_t targetLen,
                        int32_t maxDiffs, int32_t bandwidth,
                        std::shared_ptr<Alignment::SESScratchSpace> ss = nullptr)
{
    return Alignment::SES2AlignBanded<Alignment::SESAlignMode::Semiglobal,
                                      Alignment::SESTrimmingMode::Disabled,
                                      Alignment::SESTracebackMode::Enabled>(
        query, queryLen, target, targetLen, maxDiffs, bandwidth, ss);
}

auto AlignNoTraceback(const char* query, size_t queryLen, const char* target, size_t targetLen,
                      int32_t maxDiffs, int32_t bandwidth,
                      std::shared_ptr<Alignment::SESScratchSpace> ss = nullptr)
{
    return Alignment::SES2AlignBanded<Alignment::SESAlignMode::Semiglobal,
                                      Alignment::SESTrimmingMode::Disabled,
                                      Alignment::SESTracebackMode::Disabled>(
        query, queryLen, target, targetLen, maxDiffs, bandwidth, ss);
}

static const int32_t MIN_DIFFS_CAP = 10;
static const int32_t MIN_BANDWIDTH_CAP = 10;
static const int32_t MASK_DEGREE = 3;

MapperResult Mapper::Map(const PacBio::Pancake::SeqDBReaderCachedBlock& targetSeqs,
                         const PacBio::Pancake::SeedIndex& index,
                         const PacBio::Pancake::FastaSequenceCached& querySeq,
                         const PacBio::Pancake::SequenceSeedsCached& querySeeds, int64_t freqCutoff,
                         bool generateFlippedOverlap) const
{
#ifdef PANCAKE_DEBUG
    PBLOG_INFO << "Mapping query ID = " << querySeq.Id() << ", header = " << querySeq.Name();
#endif

    if (querySeq.Size() < settings_.MinQueryLen) {
        return {};
    }

    TicToc ttCollectHits;
    std::vector<SeedHit> hits;
    index.CollectHits(querySeeds.Seeds(), querySeeds.Size(), hits, freqCutoff);
    ttCollectHits.Stop();

    TicToc ttSortHits;
    std::sort(hits.begin(), hits.end(), [](const auto& a, const auto& b) {
        return PackSeedHitWithDiagonalTo128_(a) < PackSeedHitWithDiagonalTo128_(b);
    });
    ttSortHits.Stop();

    // PBLOG_INFO << "Hits: " << hits.size();

    TicToc ttChain;
    auto overlaps = FormDiagonalAnchors_(hits, querySeq, index.GetCache(), settings_.ChainBandwidth,
                                         settings_.MinNumSeeds, settings_.MinChainSpan, true,
                                         settings_.SkipSymmetricOverlaps);
    ttChain.Stop();
#ifdef PANCAKE_DEBUG
    PBLOG_INFO << "Formed diagonal anchors: " << overlaps.size();
#endif

    // Filter out multiple hits per query-target pair (e.g. tandem repeats) by
    // taking only the longest overlap chain.
    TicToc ttFilterTandem;
    if (settings_.OneHitPerTarget) {
        overlaps = FilterTandemOverlaps_(overlaps);
    }
    ttFilterTandem.Stop();
#ifdef PANCAKE_DEBUG
    PBLOG_INFO << "Overlaps after tandem filtering: " << overlaps.size();
#endif

    TicToc ttAlign;
    overlaps =
        AlignOverlaps_(targetSeqs, querySeq, overlaps, settings_.AlignmentBandwidth,
                       settings_.AlignmentMaxD, settings_.UseTraceback, settings_.NoSNPsInIdentity,
                       settings_.NoIndelsInIdentity, settings_.MaskHomopolymers,
                       settings_.MaskSimpleRepeats, generateFlippedOverlap, sesScratch_);
    ttAlign.Stop();
#ifdef PANCAKE_DEBUG
    PBLOG_INFO << "Overlaps after alignment: " << overlaps.size();
    for (const auto& ovl : overlaps) {
        PBLOG_INFO << OverlapWriterBase::PrintOverlapAsM4(ovl, "", "", true, false);
    }
#endif

    TicToc ttFilter;
    overlaps = FilterOverlaps_(overlaps, settings_.MinNumSeeds, settings_.MinIdentity,
                               settings_.MinMappedLength, settings_.MinQueryLen,
                               settings_.MinTargetLen, settings_.AllowedDovetailDist,
                               settings_.AllowedHeuristicExtendDist, settings_.BestN);
    ttFilter.Stop();

#ifdef PANCAKE_DEBUG
    PBLOG_INFO << "Overlaps after filtering: " << overlaps.size();
    for (const auto& ovl : overlaps) {
        OverlapWriterBase::PrintOverlapAsM4(stderr, ovl, querySeq.Name(),
                                            targetSeqs.GetSequence(ovl->Bid).Name(), false, false);
    }
    PBLOG_INFO << "Num anchors: " << overlaps.size();
    PBLOG_INFO << "Collected " << hits.size() << " hits.";
    PBLOG_INFO << "Time - collecting hits: " << ttCollectHits.GetMillisecs() << " ms / "
               << ttCollectHits.GetCpuMillisecs() << " CPU ms";
    PBLOG_INFO << "Time - sorting: " << ttSortHits.GetMillisecs() << " ms / "
               << ttSortHits.GetCpuMillisecs() << " CPU ms";
    PBLOG_INFO << "Time - chaining: " << ttChain.GetMillisecs() << " ms / "
               << ttChain.GetCpuMillisecs() << " CPU ms";
    PBLOG_INFO << "Time - tandem filter: " << ttFilterTandem.GetMillisecs() << " ms / "
               << ttFilterTandem.GetCpuMillisecs() << " CPU ms";
    PBLOG_INFO << "Time - alignment: " << ttAlign.GetMillisecs() << " ms / "
               << ttAlign.GetCpuMillisecs() << " CPU ms";
    PBLOG_INFO << "Time - filter: " << ttFilter.GetMillisecs() << " ms / "
               << ttFilter.GetCpuMillisecs() << " CPU ms";
    DebugWriteSeedHits_("temp/debug/mapper-0-seed_hits.csv", hits, 30, querySeq.Name(),
                        querySeq.Size(), "target", 0);
#endif

    MapperResult result;
    std::swap(result.overlaps, overlaps);
    return result;
}

OverlapPtr Mapper::MakeOverlap_(const std::vector<SeedHit>& sortedHits,
                                const PacBio::Pancake::FastaSequenceCached& querySeq,
                                const std::shared_ptr<PacBio::Pancake::SeedDBIndexCache> indexCache,
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

    const auto& sl = indexCache->GetSeedsLine(targetId);
    const int32_t targetLen = sl.numBases;

    return createOverlap(querySeq.Id(), targetId, score, identity, 0, beginHit.queryPos,
                         endHit.queryPos, querySeq.Size(), beginHit.targetRev, beginHit.targetPos,
                         endHit.targetPos, targetLen, editDist, numSeeds, OverlapType::Unknown,
                         OverlapType::Unknown);
}

std::vector<OverlapPtr> Mapper::FormDiagonalAnchors_(
    const std::vector<SeedHit>& sortedHits, const PacBio::Pancake::FastaSequenceCached& querySeq,
    const std::shared_ptr<PacBio::Pancake::SeedDBIndexCache> indexCache, int32_t chainBandwidth,
    int32_t minNumSeeds, int32_t minChainSpan, bool skipSelfHits, bool skipSymmetricOverlaps)
{

    if (sortedHits.empty()) {
        return {};
    }

    std::vector<OverlapPtr> overlaps;

    const int32_t numHits = static_cast<int32_t>(sortedHits.size());
    int32_t beginId = 0;
    int32_t beginDiag = sortedHits[beginId].Diagonal();

    // This is a combination of <targetPos, queryPos>, intended for simple comparison
    // without defining a custom complicated comparison operator.
    uint64_t minTargetQueryPosCombo = (static_cast<uint64_t>(sortedHits[beginId].targetPos) << 32) |
                                      (static_cast<uint64_t>(sortedHits[beginId].queryPos));
    uint64_t maxTargetQueryPosCombo = minTargetQueryPosCombo;
    int32_t minPosId = 0;
    int32_t maxPosId = 0;

    for (int32_t i = 0; i < numHits; ++i) {
        const auto& prevHit = sortedHits[beginId];
        const auto& currHit = sortedHits[i];
        const int32_t currDiag = currHit.Diagonal();
        const int32_t diagDiff = abs(currDiag - beginDiag);
        const uint64_t targetQueryPosCombo =
            (static_cast<uint64_t>(sortedHits[i].targetPos) << 32) |
            (static_cast<uint64_t>(sortedHits[i].queryPos));

#ifdef PANCAKE_DEBUG
        PBLOG_INFO << "[hit " << i << "] tid = " << currHit.targetId
                   << ", trev = " << currHit.targetRev << ", tpos = " << currHit.targetPos
                   << ", qpos = " << currHit.queryPos << ", flag = " << currHit.flags
                   << "; minPosId = " << minPosId << ", maxPosId = " << maxPosId;
#endif

        if (currHit.targetId != prevHit.targetId || currHit.targetRev != prevHit.targetRev ||
            diagDiff > chainBandwidth) {
            auto ovl =
                MakeOverlap_(sortedHits, querySeq, indexCache, beginId, i, minPosId, maxPosId);
            beginId = i;
            beginDiag = currDiag;

#ifdef PANCAKE_DEBUG
            PBLOG_INFO << "ovl->NumSeeds = " << ovl->NumSeeds << " (" << minNumSeeds
                       << "), minChainSpan = " << minChainSpan
                       << ", ovl->ASpan() = " << ovl->ASpan() << ", ovl->BSpan() = " << ovl->BSpan()
                       << ", skipSelfHits = " << skipSelfHits << ", ovl->Aid = " << ovl->Aid
                       << ", ovl->Bid = " << ovl->Bid
                       << ", skipSymmetricOverlaps = " << skipSymmetricOverlaps;
            PBLOG_INFO << OverlapWriterBase::PrintOverlapAsM4(ovl, "", "", true, false);
#endif

            // Add a new overlap.
            if (ovl->NumSeeds >= minNumSeeds && ovl->ASpan() > minChainSpan &&
                ovl->BSpan() > minChainSpan &&
                (skipSelfHits == false || (skipSelfHits && ovl->Bid != ovl->Aid)) &&
                (skipSymmetricOverlaps == false ||
                 (skipSymmetricOverlaps && ovl->Bid < ovl->Aid))) {
                overlaps.emplace_back(std::move(ovl));
            }
            minPosId = maxPosId = i;
            minTargetQueryPosCombo = maxTargetQueryPosCombo = targetQueryPosCombo;
        }

        // Track the minimum and maximum target positions for each diagonal.
        if (targetQueryPosCombo < minTargetQueryPosCombo) {
            minPosId = i;
            minTargetQueryPosCombo = targetQueryPosCombo;
        }
        if (targetQueryPosCombo > maxTargetQueryPosCombo) {
            maxPosId = i;
            maxTargetQueryPosCombo = targetQueryPosCombo;
        }
    }

    if ((numHits - beginId) > 0) {
        auto ovl =
            MakeOverlap_(sortedHits, querySeq, indexCache, beginId, numHits, minPosId, maxPosId);

#ifdef PANCAKE_DEBUG
        PBLOG_INFO << "ovl->NumSeeds = " << ovl->NumSeeds << " (" << minNumSeeds
                   << "), minChainSpan = " << minChainSpan << ", ovl->ASpan() = " << ovl->ASpan()
                   << ", ovl->BSpan() = " << ovl->BSpan() << ", skipSelfHits = " << skipSelfHits
                   << ", ovl->Aid = " << ovl->Aid << ", ovl->Bid = " << ovl->Bid
                   << ", skipSymmetricOverlaps = " << skipSymmetricOverlaps;
#endif

        // Add a new overlap.
        if (ovl->NumSeeds >= minNumSeeds && ovl->ASpan() > minChainSpan &&
            ovl->BSpan() > minChainSpan &&
            (skipSelfHits == false || (skipSelfHits && ovl->Bid != ovl->Aid)) &&
            (skipSymmetricOverlaps == false || (skipSymmetricOverlaps && ovl->Bid < ovl->Aid))) {

            overlaps.emplace_back(std::move(ovl));
        }
    }

    return overlaps;
}

std::vector<OverlapPtr> Mapper::FilterOverlaps_(const std::vector<OverlapPtr>& overlaps,
                                                int32_t minNumSeeds, float minIdentity,
                                                int32_t minMappedSpan, int32_t minQueryLen,
                                                int32_t minTargetLen, int32_t allowedDovetailDist,
                                                int32_t allowedExtendDist, int32_t bestN)
{

    std::vector<OverlapPtr> ret;

    for (const auto& ovl : overlaps) {
        if (100 * ovl->Identity < minIdentity || ovl->ASpan() < minMappedSpan ||
            ovl->BSpan() < minMappedSpan || ovl->NumSeeds < minNumSeeds ||
            ovl->Alen < minQueryLen || ovl->Blen < minTargetLen) {
            continue;
        }
        auto newOvl = createOverlap(ovl);
        newOvl->Atype =
            DetermineOverlapType(ovl->Arev, ovl->AstartFwd(), ovl->AendFwd(), ovl->Alen, ovl->Brev,
                                 ovl->BstartFwd(), ovl->BendFwd(), ovl->Blen, allowedDovetailDist);
        // Arev and Brev are intentionally out of place here! The A-read's orientation should always
        // be FWD, so the B-read is the one that determines the direction.
        newOvl->Btype =
            DetermineOverlapType(ovl->Arev, ovl->BstartFwd(), ovl->BendFwd(), ovl->Blen, ovl->Brev,
                                 ovl->AstartFwd(), ovl->AendFwd(), ovl->Alen, allowedDovetailDist);
        HeuristicExtendOverlapFlanks(newOvl, allowedExtendDist);
        ret.emplace_back(std::move(newOvl));
    }

    // Score is negative, as per legacy Falcon convention.
    std::stable_sort(ret.begin(), ret.end(),
                     [](const auto& a, const auto& b) { return a->Score < b->Score; });

    int32_t nToKeep = (bestN > 0) ? std::min(bestN, static_cast<int32_t>(ret.size())) : ret.size();
    ret.resize(nToKeep);

    return ret;
}

std::vector<OverlapPtr> Mapper::FilterTandemOverlaps_(const std::vector<OverlapPtr>& overlaps)
{
    if (overlaps.empty()) {
        return {};
    }

    // Make an internal copy for sorting.
    std::vector<OverlapPtr> overlaps_copy;
    for (const auto& ovl : overlaps) {
        overlaps_copy.emplace_back(createOverlap(ovl));
    }

    // Sort by length.
    std::sort(overlaps_copy.begin(), overlaps_copy.end(), [](const auto& a, const auto& b) {
        return a->Bid < b->Bid ||
               (a->Bid == b->Bid &&
                std::max(a->ASpan(), a->BSpan()) > std::max(b->ASpan(), b->BSpan()));
    });

    // Accumulate the results.
    std::vector<OverlapPtr> ret;
    ret.emplace_back(createOverlap(overlaps_copy.front()));
    for (const auto& ovl : overlaps_copy) {
        if (ovl->Bid == ret.back()->Bid) {
            continue;
        }
        ret.emplace_back(createOverlap(ovl));
    }

    return ret;
}

std::vector<OverlapPtr> Mapper::AlignOverlaps_(
    const PacBio::Pancake::SeqDBReaderCachedBlock& targetSeqs,
    const PacBio::Pancake::FastaSequenceCached& querySeq, const std::vector<OverlapPtr>& overlaps,
    double alignBandwidth, double alignMaxDiff, bool useTraceback, bool noSNPs, bool noIndels,
    bool maskHomopolymers, bool maskSimpleRepeats, bool generateFlippedOverlap,
    std::shared_ptr<PacBio::Pancake::Alignment::SESScratchSpace> sesScratch)
{
    std::vector<OverlapPtr> ret;
    const std::string reverseQuerySeq =
        PacBio::Pancake::ReverseComplement(querySeq.Bases(), querySeq.Size(), 0, querySeq.Size());

    for (size_t i = 0; i < overlaps.size(); ++i) {
        const auto& targetSeq = targetSeqs.GetSequence(overlaps[i]->Bid);
        OverlapPtr newOverlap = AlignOverlap_(
            targetSeq, querySeq, reverseQuerySeq, overlaps[i], alignBandwidth, alignMaxDiff,
            useTraceback, noSNPs, noIndels, maskHomopolymers, maskSimpleRepeats, sesScratch);

        OverlapPtr newOverlapFlipped = nullptr;
        if (generateFlippedOverlap) {
            newOverlapFlipped =
                GenerateFlippedOverlap_(targetSeq, querySeq, reverseQuerySeq, newOverlap, noSNPs,
                                        noIndels, maskHomopolymers, maskSimpleRepeats);
        }

        if (generateFlippedOverlap && newOverlap != nullptr && newOverlapFlipped == nullptr) {
            throw std::runtime_error(
                "Problem generating the flipped overlap, it's nullptr! Before flipping: " +
                OverlapWriterBase::PrintOverlapAsM4(newOverlap, "", "", true, false));
        }

        if (newOverlap != nullptr) {
            ret.emplace_back(std::move(newOverlap));
        }
        if (newOverlapFlipped != nullptr) {
            ret.emplace_back(std::move(newOverlapFlipped));
        }
    }

    return ret;
}

std::string Mapper::FetchTargetSubsequence_(const PacBio::Pancake::FastaSequenceCached& targetSeq,
                                            int32_t seqStart, int32_t seqEnd, bool revCmp)
{
    return FetchTargetSubsequence_(targetSeq.Bases(), targetSeq.Size(), seqStart, seqEnd, revCmp);
}

std::string Mapper::FetchTargetSubsequence_(const char* seq, int32_t seqLen, int32_t seqStart,
                                            int32_t seqEnd, bool revCmp)
{
    if (seqEnd == seqStart) {
        return {};
    }
    if (seqStart < 0 || seqEnd < 0 || seqStart > seqLen || seqEnd > seqLen || seqEnd < seqStart) {
        std::ostringstream oss;
        oss << "Invalid seqStart or seqEnd in a call to FetchTargetSubsequence_. seqStart = "
            << seqStart << ", seqEnd = " << seqEnd << ", seqLen = " << seqLen
            << ", revCmp = " << revCmp << ".";
        std::cerr << oss.str() << "\n";
        throw std::runtime_error(oss.str());
    }
    seqEnd = (seqEnd == 0) ? seqLen : seqEnd;
    std::string ret(seq + seqStart, seqEnd - seqStart);
    if (revCmp) {
        ret = PacBio::Pancake::ReverseComplement(ret, 0, ret.size());
    }
    return ret;
}

OverlapPtr Mapper::AlignOverlap_(
    const PacBio::Pancake::FastaSequenceCached& targetSeq,
    const PacBio::Pancake::FastaSequenceCached& querySeq, const std::string reverseQuerySeq,
    const OverlapPtr& ovl, double alignBandwidth, double alignMaxDiff, bool useTraceback,
    bool noSNPs, bool noIndels, bool maskHomopolymers, bool maskSimpleRepeats,
    std::shared_ptr<PacBio::Pancake::Alignment::SESScratchSpace> sesScratch)
{

    if (ovl == nullptr) {
        return nullptr;
    }

#ifdef PANCAKE_DEBUG_ALN
    PBLOG_INFO << "Initial: " << OverlapWriterBase::PrintOverlapAsM4(ovl, "", "", true, false);
#endif

    OverlapPtr ret = createOverlap(ovl);
    PacBio::Pancake::Alignment::SesResults sesResultRight;
    PacBio::Pancake::Alignment::SesResults sesResultLeft;

    ///////////////////////////
    /// Align forward pass. ///
    ///////////////////////////
    {
        const int32_t qStart = ovl->Astart;
        const int32_t qEnd = ovl->Alen;
        const int32_t qSpan = qEnd - qStart;
        const int32_t tStartFwd = ovl->Brev ? (ovl->Blen - ovl->Bend) : ovl->Bstart;
        const int32_t tEndFwd = ovl->Brev ? (ovl->Blen - ovl->Bstart) : ovl->Bend;
        std::string tseq;
        if (ovl->Brev) {
            // Extract reverse complemented target sequence.
            // The reverse complement begins at the last mapped position (tEndFwd),
            // and ends at the the tStartFwd reduced by an allowed overhang.
            // The allowed overhang is designed to absorb the entire unmapped end
            // of the query, and at most 2x that length in the target. (The 2xhang
            // is just to be safe, because no proper alignment should be 2x longer
            // in one sequence than in the other).
            int32_t minHangLen = std::min(ovl->Alen - ovl->Aend, tStartFwd);
            int32_t extractBegin = std::max(0, tStartFwd - minHangLen * 2);
            int32_t extractEnd = tEndFwd;
            tseq = FetchTargetSubsequence_(targetSeq, extractBegin, extractEnd, ovl->Brev);
        } else {
            // Take the sequence starting from the start position, and reaching
            // until the end of the query (or target, which ever is the shorter).
            // Extract 2x larger overhang to be safe - no proper alignment should be
            // 2x longer in one sequence than the other.
            int32_t minHangLen = std::min(ovl->Blen - tEndFwd, ovl->Alen - ovl->Aend);
            int32_t extractBegin = tStartFwd;
            int32_t extractEnd = std::min(ovl->Blen, tEndFwd + minHangLen * 2);
            tseq = FetchTargetSubsequence_(targetSeq, extractBegin, extractEnd, ovl->Brev);
        }
        const int32_t tSpan = tseq.size();
        const int32_t dMax =
            std::max(MIN_DIFFS_CAP, static_cast<int32_t>(ovl->Alen * alignMaxDiff));
        const int32_t bandwidth =
            std::max(MIN_BANDWIDTH_CAP,
                     static_cast<int32_t>(std::min(ovl->Blen, ovl->Alen) * alignBandwidth));

        if (useTraceback) {
            sesResultRight = AlignWithTraceback(querySeq.Bases() + qStart, qSpan, tseq.c_str(),
                                                tSpan, dMax, bandwidth, sesScratch);
        } else {
            sesResultRight = AlignNoTraceback(querySeq.Bases() + qStart, qSpan, tseq.c_str(), tSpan,
                                              dMax, bandwidth, sesScratch);
        }

        ret->Aend = sesResultRight.lastQueryPos;
        ret->Bend = sesResultRight.lastTargetPos;
        ret->Aend += ovl->Astart;
        ret->Bend += ovl->Bstart;
        ret->EditDistance = sesResultRight.numDiffs;
        ret->Score = -(std::min(ret->ASpan(), ret->BSpan()) - ret->EditDistance);

#ifdef PANCAKE_DEBUG_ALN
        PBLOG_INFO << "dMax = " << dMax << ", bandwidth = " << bandwidth;
        PBLOG_INFO << "Right: diffs = " << sesResultRight.diffs;
        PBLOG_INFO << "After right: "
                   << OverlapWriterBase::PrintOverlapAsM4(ret, "", "", true, false);
#endif
    }

    ///////////////////////////
    /// Align reverse pass. ///
    ///////////////////////////
    {
        // Reverse query coordinates.
        const int32_t qStart = ret->Alen - ret->Astart;
        const int32_t qEnd = ret->Alen;
        const int32_t qSpan = qEnd - qStart;
        const int32_t tStartFwd = ret->Brev ? (ret->Blen - ret->Bend) : ret->Bstart;
        const int32_t tEndFwd = ret->Brev ? (ret->Blen - ret->Bstart) : ret->Bend;
        std::string tseq;
        if (ovl->Brev) {
            int32_t minHangLen = std::min(ovl->Blen - tEndFwd, qStart);
            int32_t extractBegin = tEndFwd;
            int32_t extractEnd = std::min(ret->Blen, tEndFwd + minHangLen * 2);
            tseq = FetchTargetSubsequence_(targetSeq, extractBegin, extractEnd, !ret->Brev);
        } else {
            int32_t minHangLen = std::min(ovl->Astart, tStartFwd);
            int32_t extractBegin = std::max(0, tStartFwd - minHangLen * 2);
            int32_t extractEnd = tStartFwd;
            tseq = FetchTargetSubsequence_(targetSeq, extractBegin, extractEnd, !ret->Brev);
        }
        const int32_t tSpan = tseq.size();
        const int32_t dMax = std::max(MIN_DIFFS_CAP, static_cast<int32_t>(ovl->Alen * alignMaxDiff -
                                                                          sesResultRight.numDiffs));
        const int32_t bandwidth =
            std::max(MIN_BANDWIDTH_CAP,
                     static_cast<int32_t>(std::min(ovl->Blen, ovl->Alen) * alignBandwidth));

        if (useTraceback) {
            sesResultLeft = AlignWithTraceback(reverseQuerySeq.c_str() + qStart, qSpan,
                                               tseq.c_str(), tSpan, dMax, bandwidth, sesScratch);
        } else {
            sesResultLeft = AlignNoTraceback(reverseQuerySeq.c_str() + qStart, qSpan, tseq.c_str(),
                                             tSpan, dMax, bandwidth, sesScratch);
        }

        ret->Astart = ovl->Astart - sesResultLeft.lastQueryPos;
        ret->Bstart = ovl->Bstart - sesResultLeft.lastTargetPos;
        std::reverse(sesResultLeft.cigar.begin(), sesResultLeft.cigar.end());
    }

    PacBio::Pancake::Alignment::DiffCounts diffs =
        sesResultRight.diffCounts + sesResultLeft.diffCounts;

    // Compute edit distance, identity and score.
    ret->EditDistance = -1;
    if (useTraceback) {
        diffs.Identity(noSNPs, noIndels, ret->Identity, ret->EditDistance);
        ret->Score = -diffs.numEq;
    } else {
        ret->EditDistance = sesResultLeft.numDiffs + sesResultRight.numDiffs;
        ret->Score = -(std::min(ret->ASpan(), ret->BSpan()) - ret->EditDistance);
        const float span = std::max(ret->ASpan(), ret->BSpan());
        ret->Identity =
            ((span != 0.0f) ? ((span - static_cast<float>(ret->EditDistance)) / span) : -0.0f);
    }

    // Merge the CIGAR strings.
    ret->Cigar = std::move(sesResultLeft.cigar);
    if (sesResultRight.cigar.size() > 0) {
        AppendToCigar(ret->Cigar, sesResultRight.cigar.front().Type(),
                      sesResultRight.cigar.front().Length());
        ret->Cigar.insert(ret->Cigar.end(), sesResultRight.cigar.begin() + 1,
                          sesResultRight.cigar.end());
    }

    NormalizeAndExtractVariantsInPlace_(ret, targetSeq, querySeq, reverseQuerySeq, noSNPs, noIndels,
                                        maskHomopolymers, maskSimpleRepeats);

#ifdef PANCAKE_DEBUG_ALN
    PBLOG_INFO << "Final: " << OverlapWriterBase::PrintOverlapAsM4(ret, "", "", true, false);
#endif

    return ret;
}

OverlapPtr Mapper::GenerateFlippedOverlap_(const PacBio::Pancake::FastaSequenceCached& targetSeq,
                                           const PacBio::Pancake::FastaSequenceCached& querySeq,
                                           const std::string reverseQuerySeq, const OverlapPtr& ovl,
                                           bool noSNPs, bool noIndels, bool maskHomopolymers,
                                           bool maskSimpleRepeats)
{

    if (ovl == nullptr) {
        return nullptr;
    }

    auto ret = CreateFlippedOverlap(ovl);

    NormalizeAndExtractVariantsInPlace_(ret, targetSeq, querySeq, reverseQuerySeq, noSNPs, noIndels,
                                        maskHomopolymers, maskSimpleRepeats);

    return ret;
}

void Mapper::NormalizeAndExtractVariantsInPlace_(
    OverlapPtr& ovl, const PacBio::Pancake::FastaSequenceCached& targetSeq,
    const PacBio::Pancake::FastaSequenceCached& querySeq, const std::string reverseQuerySeq,
    bool noSNPs, bool noIndels, bool maskHomopolymers, bool maskSimpleRepeats)
{
    // Extract the variant strings.
    if (ovl->Cigar.empty()) {
        return;
    }

    const char* Aseq = querySeq.Bases();
    int32_t Alen = querySeq.Size();
    const char* Bseq = targetSeq.Bases();
    int32_t Blen = targetSeq.Size();

    if (ovl->IsFlipped) {
        Aseq = targetSeq.Bases();
        Alen = targetSeq.Size();
        Bseq = querySeq.Bases();
        Blen = querySeq.Size();
    }

    // / This works, but is slower because it requires to copy the target sequence every time.
    // / Since we already have the reversed query, we can simply reverse the CIGAR string and provide
    // / the reversed query, as below.
    PacBio::Pancake::Alignment::DiffCounts diffsPerBase;
    PacBio::Pancake::Alignment::DiffCounts diffsPerEvent;

    auto tseq = FetchTargetSubsequence_(Bseq, Blen, ovl->BstartFwd(), ovl->BendFwd(), ovl->Brev);

    const char* querySub = Aseq + ovl->Astart;
    int32_t querySubLen = ovl->ASpan();
    const char* targetSub = tseq.c_str();
    int32_t targetSubLen = tseq.size();
    auto& cigar = ovl->Cigar;

    // If the B-read is reversed, then we need to reverse and left-align the CIGAR string
    // to be consistent with the other alignments which were in forward.
    // The reason is that in the rest of the workflow up to this point we treated the
    // query as either fwd or reverse, but finally we actually take query as FWD in the output.
    // Normalization needs to happen in this context, so that all left-aligned indels line up.
    if (ovl->Brev) {
        auto tempCigar = ovl->Cigar;
        std::reverse(tempCigar.begin(), tempCigar.end());
        auto qseq =
            FetchTargetSubsequence_(Aseq, Alen, ovl->AstartFwd(), ovl->AendFwd(), ovl->Brev);

        tempCigar = NormalizeCigar(qseq.c_str(), ovl->ASpan(), Bseq + ovl->BstartFwd(),
                                   ovl->BSpan(), tempCigar);
        std::reverse(tempCigar.begin(), tempCigar.end());
        std::swap(cigar, tempCigar);

    } else {
        cigar = NormalizeCigar(querySub, querySubLen, targetSub, targetSubLen, cigar);
    }

    ExtractVariantString(querySub, querySubLen, targetSub, targetSubLen, cigar, maskHomopolymers,
                         maskSimpleRepeats, ovl->Avars, ovl->Bvars, diffsPerBase, diffsPerEvent);
}

void NormalizeAndExtractVariantsInPlaceDeprecated_(
    OverlapPtr& ovl, const PacBio::Pancake::FastaSequenceCached& targetSeq,
    const PacBio::Pancake::FastaSequenceCached& querySeq, const std::string reverseQuerySeq,
    bool noSNPs, bool noIndels, bool maskHomopolymers, bool maskSimpleRepeats)
{
    // Extract the variant strings.
    if (ovl->Cigar.empty()) {
        return;
    }

    /// This works, but is slower because it requires to copy the target sequence every time.
    /// Since we already have the reversed query, we can simply reverse the CIGAR string and provide
    /// the reversed query, as below.
    // auto tseq = FetchTargetSubsequence_(targetSeq, ovl->BstartFwd(), ovl->BendFwd(), ovl->Brev);
    // ovl->VariantString = ExtractVariantString(querySeq.Bases() + ovl->Astart, ovl->ASpan(), tseq.c_str(), tseq.size(),
    //               ovl->Cigar, false, false);

    PacBio::Pancake::Alignment::DiffCounts diffsPerBase;
    PacBio::Pancake::Alignment::DiffCounts diffsPerEvent;

    if (ovl->Brev) {
        // Get the correct subsequences as C-strings.
        const char* querySub = reverseQuerySeq.c_str() + (ovl->Alen - ovl->Aend);
        int32_t querySubLen = ovl->ASpan();
        const char* targetSub = targetSeq.Bases() + ovl->BstartFwd();
        int32_t targetSubLen = ovl->BSpan();
        // Reverse the CIGAR (it was generated earlier by reversing the target).
        auto tempCigar = ovl->Cigar;
        std::reverse(tempCigar.begin(), tempCigar.end());

        // Normalize the CIGAR operations.
        tempCigar = NormalizeCigar(querySub, querySubLen, targetSub, targetSubLen, tempCigar);
        ovl->Cigar = tempCigar;
        std::reverse(ovl->Cigar.begin(), ovl->Cigar.end());

        // Extract variants for the updated CIGAR.
        ExtractVariantString(querySub, querySubLen, targetSub, targetSubLen, tempCigar,
                             maskHomopolymers, maskSimpleRepeats, ovl->Avars, ovl->Bvars,
                             diffsPerBase, diffsPerEvent);

        ovl->Avars = Pancake::ReverseComplement(ovl->Avars, 0, ovl->Avars.size());
        ovl->Bvars = Pancake::ReverseComplement(ovl->Bvars, 0, ovl->Bvars.size());
    } else {
        // Get the correct subsequences as C-strings.
        const char* querySub = querySeq.Bases() + ovl->Astart;
        int32_t querySubLen = ovl->ASpan();
        const char* targetSub = targetSeq.Bases() + ovl->BstartFwd();
        int32_t targetSubLen = ovl->BSpan();
        auto& tempCigar = ovl->Cigar;
        // Normalize the CIGAR operations.
        tempCigar = NormalizeCigar(querySub, querySubLen, targetSub, targetSubLen, tempCigar);
        // Extract variants for the updated CIGAR.
        ExtractVariantString(querySub, querySubLen, targetSub, targetSubLen, tempCigar,
                             maskHomopolymers, maskSimpleRepeats, ovl->Avars, ovl->Bvars,
                             diffsPerBase, diffsPerEvent);
    }

    const auto& diffs = diffsPerBase;
    diffs.Identity(noSNPs, noIndels, ovl->Identity, ovl->EditDistance);
    ovl->Score = -diffs.numEq;
}

void Mapper::DebugWriteSeedHits_(const std::string& outPath, const std::vector<SeedHit>& hits,
                                 int32_t seedLen, const std::string& queryName, int64_t queryLen,
                                 const std::string& targetName, int64_t targetLen)
{
    std::ofstream ofs(outPath);
    // Simply walk away if the file cannot be open.
    // Avoid writing debug output if a specific path is not available.
    if (ofs.is_open() == false) {
        return;
    }
    ofs << queryName << "\t0\t" << queryLen << "\t" << targetName << "\t0\t" << targetLen << "\t0.0"
        << std::endl;
    for (size_t i = 0; i < hits.size(); ++i) {
        int32_t clusterId = hits[i].targetId * 2 + (hits[i].targetRev ? 1 : 0);
        ofs << hits[i].queryPos << "\t" << hits[i].targetPos << "\t"
            << "\t" << clusterId << std::endl;
        ofs << hits[i].queryPos + seedLen << "\t" << hits[i].targetPos + seedLen << "\t"
            << clusterId << std::endl;
    }
}

PacBio::Pancake::Int128t Mapper::PackSeedHitWithDiagonalTo128_(const SeedHit& sh)
{
    PacBio::Pancake::Int128t ret = 0;
    const int32_t diag = (sh.targetPos - sh.queryPos);
    ret = ((static_cast<PacBio::Pancake::Int128t>(sh.targetId) & MASK128_LOW32bit) << 97) |
          ((static_cast<PacBio::Pancake::Int128t>(sh.targetRev) & MASK128_LOW32bit) << 96) |
          ((static_cast<PacBio::Pancake::Int128t>(diag) & MASK128_LOW32bit) << 64) |
          ((static_cast<PacBio::Pancake::Int128t>(sh.targetPos) & MASK128_LOW32bit) << 32) |
          (static_cast<PacBio::Pancake::Int128t>(sh.queryPos) & MASK128_LOW32bit);
    return ret;
}

}  // namespace Pancake
}  // namespace PacBio
