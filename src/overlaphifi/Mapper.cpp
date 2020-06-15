// Authors: Ivan Sovic

#include <lib/kxsort/kxsort.h>
#include <pacbio/alignment/AlignmentTools.h>
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
#include <pacbio/alignment/SesAlignBanded.hpp>
#include <sstream>

namespace PacBio {
namespace Pancake {

// #define PANCAKE_DEBUG

MapperResult Mapper::Map(const PacBio::Pancake::SeqDBReaderCachedBlock& targetSeqs,
                         const PacBio::Pancake::SeedIndex& index,
                         const PacBio::Pancake::FastaSequenceCached& querySeq,
                         const PacBio::Pancake::SequenceSeedsCached& querySeeds,
                         int64_t freqCutoff) const
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

    TicToc ttChain;
    auto overlaps = FormDiagonalAnchors_(hits, querySeq, index.GetCache(), settings_.ChainBandwidth,
                                         settings_.MinNumSeeds, settings_.MinChainSpan, true,
                                         settings_.SkipSymmetricOverlaps);
    ttChain.Stop();

    // Filter out multiple hits per query-target pair (e.g. tandem repeats) by
    // taking only the longest overlap chain.
    TicToc ttFilterTandem;
    if (settings_.OneHitPerTarget) {
        overlaps = FilterTandemOverlaps_(overlaps);
    }
    ttFilterTandem.Stop();

    TicToc ttAlign;
    overlaps =
        AlignOverlaps_(targetSeqs, querySeq, overlaps, settings_.AlignmentBandwidth,
                       settings_.AlignmentMaxD, settings_.UseTraceback, settings_.NoSNPsInIdentity,
                       settings_.NoIndelsInIdentity, sesScratch_);
    ttAlign.Stop();

    TicToc ttFilter;
    overlaps = FilterOverlaps_(overlaps, settings_.MinNumSeeds, settings_.MinIdentity,
                               settings_.MinMappedLength, settings_.MinQueryLen,
                               settings_.MinTargetLen, settings_.AllowedDovetailDist,
                               settings_.AllowedHeuristicExtendDist, settings_.BestN);
    ttFilter.Stop();

#ifdef PANCAKE_DEBUG
    for (const auto& ovl : overlaps) {
        OverlapWriterBase::PrintOverlapAsM4(stdout, ovl, querySeq.Name(),
                                            targetSeqs.GetSequence(ovl->Bid).Name(), false);
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
                         endHit.targetPos, targetLen, editDist, numSeeds, OverlapType::Unknown);
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

        if (currHit.targetId != prevHit.targetId || currHit.targetRev != prevHit.targetRev ||
            diagDiff > chainBandwidth) {
            auto ovl =
                MakeOverlap_(sortedHits, querySeq, indexCache, beginId, i, minPosId, maxPosId);
            beginId = i;
            beginDiag = currDiag;

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
        newOvl->Type = DetermineOverlapType(*ovl, allowedDovetailDist);
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
    std::shared_ptr<PacBio::Pancake::Alignment::SESScratchSpace> sesScratch)
{
    std::vector<OverlapPtr> ret;
    const std::string reverseQuerySeq =
        PacBio::Pancake::ReverseComplement(querySeq.Bases(), querySeq.Size(), 0, querySeq.Size());

    for (size_t i = 0; i < overlaps.size(); ++i) {
        const auto& targetSeq = targetSeqs.GetSequence(overlaps[i]->Bid);
        auto newOverlap =
            AlignOverlap_(targetSeq, querySeq, reverseQuerySeq, overlaps[i], alignBandwidth,
                          alignMaxDiff, useTraceback, noSNPs, noIndels, sesScratch);
        if (newOverlap == nullptr) {
            continue;
        }
        ret.emplace_back(std::move(newOverlap));
    }
    return ret;
}

std::string Mapper::FetchTargetSubsequence_(const PacBio::Pancake::FastaSequenceCached& targetSeq,
                                            int32_t seqStart, int32_t seqEnd, bool revCmp)
{
    const int32_t seqLen = targetSeq.Size();
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
    std::string ret(targetSeq.Bases() + seqStart, seqEnd - seqStart);
    if (revCmp) {
        ret = PacBio::Pancake::ReverseComplement(ret, 0, ret.size());
    }
    return ret;
}

OverlapPtr Mapper::AlignOverlap_(
    const PacBio::Pancake::FastaSequenceCached& targetSeq,
    const PacBio::Pancake::FastaSequenceCached& querySeq, const std::string reverseQuerySeq,
    const OverlapPtr& ovl, double alignBandwidth, double alignMaxDiff, bool useTraceback,
    bool noSNPs, bool noIndels,
    std::shared_ptr<PacBio::Pancake::Alignment::SESScratchSpace> sesScratch)
{

    if (ovl == nullptr) {
        return nullptr;
    }

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
        const int32_t dMax = ovl->Alen * alignMaxDiff;
        const int32_t bandwidth = std::min(ovl->Blen, ovl->Alen) * alignBandwidth;

        if (useTraceback) {
            sesResultRight =
                PacBio::Pancake::Alignment::SESAlignBanded<Alignment::SESAlignMode::Semiglobal,
                                                           Alignment::SESTracebackMode::Enabled>(
                    querySeq.Bases() + qStart, qSpan, tseq.c_str(), tSpan, dMax, bandwidth,
                    sesScratch);
        } else {
            sesResultRight = PacBio::Pancake::Alignment::SESDistanceBanded(
                querySeq.Bases() + qStart, qSpan, tseq.c_str(), tSpan, dMax, bandwidth);
        }

        ret->Aend = sesResultRight.lastQueryPos;
        ret->Bend = sesResultRight.lastTargetPos;
        ret->Aend += ovl->Astart;
        ret->Bend += ovl->Bstart;
        ret->EditDistance = sesResultRight.numDiffs;
        ret->Score = -(std::min(ret->ASpan(), ret->BSpan()) - ret->EditDistance);
        // std::cerr << "CIGAR right: " << cigarRight.ToStdString() << "\n";
        // std::reverse
        // std::cerr << "(fwd) sesResult: " << sesResult << "\n";
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
        const int32_t dMax = ovl->Alen * alignMaxDiff - sesResultRight.numDiffs;
        const int32_t bandwidth = std::min(ovl->Blen, ovl->Alen) * alignBandwidth;

        if (useTraceback) {
            sesResultLeft =
                PacBio::Pancake::Alignment::SESAlignBanded<Alignment::SESAlignMode::Semiglobal,
                                                           Alignment::SESTracebackMode::Enabled>(
                    reverseQuerySeq.c_str() + qStart, qSpan, tseq.c_str(), tSpan, dMax, bandwidth,
                    sesScratch);
        } else {
            sesResultLeft = PacBio::Pancake::Alignment::SESDistanceBanded(
                reverseQuerySeq.c_str() + qStart, qSpan, tseq.c_str(), tSpan, dMax, bandwidth);
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

    ret->Cigar = std::move(sesResultLeft.cigar);
    if (sesResultRight.cigar.size() > 0) {
        AppendToCigar(ret->Cigar, sesResultRight.cigar.front().Type(),
                      sesResultRight.cigar.front().Length());
        ret->Cigar.insert(ret->Cigar.end(), sesResultRight.cigar.begin() + 1,
                          sesResultRight.cigar.end());
    }

    return ret;
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
