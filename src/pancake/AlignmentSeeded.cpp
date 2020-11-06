// Authors: Ivan Sovic

#include <pacbio/alignment/AlignmentTools.h>
#include <pacbio/pancake/AlignmentSeeded.h>
#include <pbcopper/logging/Logging.h>

namespace PacBio {
namespace Pancake {

// #define DEBUG_ALIGNMENT_SEEDED

RegionsToAlign ExtractAlignmentRegions(const std::vector<SeedHit>& inSortedHits,
                                       const char* querySeqFwd, const char* querySeqRev,
                                       int32_t qLen, const char* targetSeq, int32_t tLen,
                                       bool isRev, int32_t minAlignmentSpan,
                                       int32_t maxFlankExtensionDist)
{
    if (inSortedHits.empty()) {
        return {};
    }

    // NOTE: This is only required if the hit coordinates are always
    // fwd in the query sequence. If they are not, this can be removed.
    //
    // Reverse the hit coordinates if required, to make the
    // alignment simpler. This is required because the original coordinate
    // system keeps query in fwd and target in strand, whereas here we
    // use query in strand and target in fwd.
    std::vector<SeedHit> hits = inSortedHits;
    if (isRev) {
        for (size_t i = 0; i < hits.size(); ++i) {
            auto& hit = hits[i];
            // std::swap(hit.queryPos, hit.targetPos
            hit.queryPos = qLen - hit.queryPos;
            hit.targetPos = tLen - hit.targetPos;

            if (i > 0) {
                bool isLongCurr = hits[i].CheckFlagLongJoin();
                bool isLongPrev = hits[i - 1].CheckFlagLongJoin();
                hits[i - 1].SetFlagLongJoin(isLongCurr);
                hits[i].SetFlagLongJoin(isLongPrev);
            }
        }
        std::reverse(hits.begin(), hits.end());
    }

    const int32_t strandId = isRev;
    const char* querySeqInStrand = (strandId == 0) ? querySeqFwd : querySeqRev;
    const char* targetSeqInStrand = targetSeq;
    const double flankExtFactor = 1.3;

    RegionsToAlign ret;

    // Align the front.
    ret.actualQueryStart = hits.front().queryPos;
    ret.actualTargetStart = hits.front().targetPos;
    if (ret.actualQueryStart > 0 && ret.actualTargetStart > 0) {
        // Determine the maximum flank we want to align, in both target and query.
        const int32_t projTStart =
            std::max(0.0, ret.actualTargetStart - flankExtFactor * ret.actualQueryStart);
        const int32_t projQStart = 0;
        const int32_t qExtLen = std::min(ret.actualQueryStart - projQStart, maxFlankExtensionDist);
        const int32_t tExtLen = std::min(ret.actualTargetStart - projTStart, maxFlankExtensionDist);

        // Find the frame of the sequences to compare.
        int32_t qStart = ret.actualQueryStart - qExtLen;
        int32_t tStart = ret.actualTargetStart - tExtLen;
        const int32_t qSpan = ret.actualQueryStart - qStart;
        const int32_t tSpan = ret.actualTargetStart - tStart;

        // Extract the sequence.
        ret.frontSemiglobal.qSeq = std::string(querySeqInStrand + qStart, qSpan);
        ret.frontSemiglobal.tSeq = std::string(targetSeqInStrand + tStart, tSpan);

        // Reverse the sequences for alignment.
        std::reverse(ret.frontSemiglobal.qSeq.begin(), ret.frontSemiglobal.qSeq.end());
        std::reverse(ret.frontSemiglobal.tSeq.begin(), ret.frontSemiglobal.tSeq.end());
    }

    // Align between seeds.
    const int32_t nHits = hits.size();
    std::vector<Data::Cigar> cigarChunks;
    int32_t startId = 0;
    for (int32_t i = 1; i < nHits; ++i) {
        // Shorthands.
        const auto& h1 = hits[startId];
        const auto& h2 = hits[i];

        // Compute the new region.
        AlignmentRegion region;
        region.qSeq = querySeqInStrand;
        region.qStart = h1.queryPos;
        region.qSpan = h2.queryPos - h1.queryPos;
        region.tSeq = targetSeq;
        region.tStart = h1.targetPos;
        region.tSpan = h2.targetPos - h1.targetPos;
        region.semiglobal = false;

        // Sanity check.
        if (region.qSpan < 0 || region.tSpan < 0) {
            std::ostringstream oss;
            oss << "qStart = " << region.qStart << ", qSpan = " << region.qSpan
                << ", tStart = " << region.tStart << ", tSpan = " << region.tSpan << "\n";
            throw std::runtime_error(oss.str());
        }

        // Skip short alignment portions.
        if ((i + 1) < nHits && h2.CheckFlagLongJoin() == false &&
            (region.qSpan < minAlignmentSpan || region.tSpan < minAlignmentSpan)) {
            continue;
        }

        // Update the start ID for the next iteration.
        startId = i;

        // Add the new region.
        ret.internalGlobal.emplace_back(std::move(region));
    }

    // Back chunk.
    ret.actualQueryEnd = hits.back().queryPos;
    ret.actualTargetEnd = hits.back().targetPos;
    if (ret.actualQueryEnd < qLen && ret.actualTargetEnd < tLen) {
        // Determine the maximum flank we want to align, in both target and query.
        const int32_t qFlankLen = qLen - ret.actualQueryEnd;
        const int32_t projTEnd =
            std::min(tLen, static_cast<int32_t>(ret.actualTargetEnd + flankExtFactor * qFlankLen));
        const int32_t projQEnd = qLen;
        const int32_t qExtLen = std::min(projQEnd - ret.actualQueryEnd, maxFlankExtensionDist);
        const int32_t tExtLen = std::min(projTEnd - ret.actualTargetEnd, maxFlankExtensionDist);

        // Compute the new region.
        AlignmentRegion region;
        region.qSeq = querySeqInStrand;
        region.qStart = ret.actualQueryEnd;
        region.qSpan = qExtLen;
        region.tSeq = targetSeqInStrand;
        region.tStart = ret.actualTargetEnd;
        region.tSpan = tExtLen;
        region.semiglobal = true;

        // Assign the region to the return object.
        ret.backSemiglobal = std::move(region);
    }

    return ret;
}

RegionsToAlignResults AlignRegionsGeneric(const RegionsToAlign& regions,
                                          AlignerBasePtr& alignerGlobal, AlignerBasePtr& alignerExt)
{
    RegionsToAlignResults ret;

    // Align the front.
    if (regions.frontSemiglobal.qSeq.size() > 0 && regions.frontSemiglobal.tSeq.size() > 0) {
        const auto& region = regions.frontSemiglobal;
        AlignmentResult alnRes = alignerExt->Extend(region.qSeq.c_str(), region.qSeq.size(),
                                                    region.tSeq.c_str(), region.tSeq.size());
        std::reverse(alnRes.cigar.begin(), alnRes.cigar.end());
        if (alnRes.cigar.size() > 0) {
            ret.cigarChunks.emplace_back(std::move(alnRes.cigar));
        }
        ret.offsetFrontQuery = alnRes.lastQueryPos;
        ret.offsetFrontTarget = alnRes.lastTargetPos;
    }

    // Alignment in between seed hits.
    for (size_t i = 0; i < regions.internalGlobal.size(); ++i) {
        const auto& region = regions.internalGlobal[i];
        AlignmentResult alnRes = alignerGlobal->Global(region.qSeq + region.qStart, region.qSpan,
                                                       region.tSeq + region.tStart, region.tSpan);
        ret.cigarChunks.emplace_back(std::move(alnRes.cigar));

#ifdef DEBUG_ALIGNMENT_SEEDED
        std::cerr << "[aln region i = " << i << " / " << regions.internalGlobal.size()
                  << "] region.qStart = " << region.qStart
                  << ", region.qEnd = " << (region.qStart + region.qSpan)
                  << ", region.qSpan = " << region.qSpan << ", region.tStart = " << region.tStart
                  << ", region.tEnd = " << (region.tStart + region.tSpan)
                  << ", region.tSpan = " << region.tSpan
                  << ", CIGAR: " << ret.cigarChunks.back().ToStdString() << "\n";

        ValidateCigar(region.qSeq + region.qStart, region.qSpan, region.tSeq + region.tStart,
                      region.tSpan, ret.cigarChunks.back(), "Chunk validation.");
#endif
    }

    // Align the back.
    if (regions.backSemiglobal.qSpan > 0 && regions.backSemiglobal.tSpan > 0) {
        const auto& region = regions.backSemiglobal;
        AlignmentResult alnRes = alignerExt->Extend(region.qSeq + region.qStart, region.qSpan,
                                                    region.tSeq + region.tStart, region.tSpan);
        if (alnRes.cigar.size() > 0) {
            ret.cigarChunks.emplace_back(std::move(alnRes.cigar));
        }
        ret.offsetBackQuery = alnRes.lastQueryPos;
        ret.offsetBackTarget = alnRes.lastTargetPos;
    }

    return ret;
}

OverlapPtr AlignOverlapSeeded(const OverlapPtr& ovl, const std::vector<SeedHit>& sortedHits,
                              const char* targetSeq, const int32_t targetLen, const char* queryFwd,
                              const char* queryRev, const int32_t queryLen,
                              int32_t minAlignmentSpan, int32_t maxFlankExtensionDist,
                              AlignerBasePtr& alignerGlobal, AlignerBasePtr& alignerExt)
{
    if (ovl->Arev) {
        throw std::runtime_error("The ovl->Arev should always be false! In Align_.");
    }
    if (sortedHits.empty()) {
        throw std::runtime_error("There are no seed hits provided to Align_ for alignment.");
    }

    // Prepare the regions for alignment.
    RegionsToAlign regions =
        ExtractAlignmentRegions(sortedHits, queryFwd, queryRev, queryLen, targetSeq, targetLen,
                                ovl->Brev, minAlignmentSpan, maxFlankExtensionDist);

    // Run the alignment.
    RegionsToAlignResults alns = AlignRegionsGeneric(regions, alignerGlobal, alignerExt);

    // Process the alignment results and make a new overlap.
    OverlapPtr ret = createOverlap(ovl);
    ret->Cigar.clear();
    ret->Astart = regions.actualQueryStart - alns.offsetFrontQuery;
    ret->Aend = regions.actualQueryEnd + alns.offsetBackQuery;
    ret->Bstart = regions.actualTargetStart - alns.offsetFrontTarget;
    ret->Bend = regions.actualTargetEnd + alns.offsetBackTarget;

    // Merge the CIGAR chunks.
    for (const auto& currCigar : alns.cigarChunks) {
        if (currCigar.empty()) {
            continue;
        }
        if (ret->Cigar.empty() || ret->Cigar.back().Type() != currCigar.front().Type()) {
            ret->Cigar.emplace_back(currCigar.front());
        } else {
            ret->Cigar.back().Length(ret->Cigar.back().Length() + currCigar.front().Length());
        }
        ret->Cigar.insert(ret->Cigar.end(), currCigar.begin() + 1, currCigar.end());
    }

    // Reverse the CIGAR and the coordinates if needed.
    if (ovl->Brev) {
        // CIGAR reversal.
        std::reverse(ret->Cigar.begin(), ret->Cigar.end());

        // Reverse the query coordinates.
        std::swap(ret->Astart, ret->Aend);
        ret->Astart = ret->Alen - ret->Astart;
        ret->Aend = ret->Alen - ret->Aend;

        // Get the forward-oriented target coordinates.
        std::swap(ret->Bstart, ret->Bend);
        ret->Bstart = ret->Blen - ret->Bstart;
        ret->Bend = ret->Blen - ret->Bend;
    }

    // Set the alignment identity and edit distance.
    Alignment::DiffCounts diffs = CigarDiffCounts(ret->Cigar);
    diffs.Identity(false, false, ret->Identity, ret->EditDistance);

    // Validation. In case an alignment was dropped.
    try {
        const int32_t qstart = ret->Astart;
        const int32_t tstart = ret->Bstart;
        const char* querySeqFwd = queryFwd;
        const std::string targetSeqForValidation =
            (ret->Brev) ? PacBio::Pancake::ReverseComplement(targetSeq, 0, targetLen) : targetSeq;
        ValidateCigar(querySeqFwd + qstart, ret->ASpan(), targetSeqForValidation.c_str() + tstart,
                      ret->BSpan(), ret->Cigar, "Full length validation.");
    } catch (std::exception& e) {
        ret = nullptr;
        PBLOG_DEBUG << "[Note: Exception when aligning!] " << e.what() << "\n";
    }

    return ret;
}

}  // namespace Pancake
}  // namespace PacBio
