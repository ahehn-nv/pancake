// Authors: Ivan Sovic

#include <pacbio/alignment/AlignmentTools.h>
#include <pacbio/pancake/AlignmentSeeded.h>
#include <pbcopper/logging/Logging.h>
#include <iostream>

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

    const double flankExtFactor = 1.3;

    RegionsToAlign ret;

    ret.targetSeq = targetSeq;
    ret.querySeqFwd = querySeqFwd;
    ret.querySeqRev = querySeqRev;
    ret.targetLen = tLen;
    ret.queryLen = qLen;

    ret.actualQueryStart = hits.front().queryPos;
    ret.actualTargetStart = hits.front().targetPos;
    ret.actualQueryEnd = hits.back().queryPos;
    ret.actualTargetEnd = hits.back().targetPos;

    int32_t numRegions = 0;

    // Extract the front.
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

        AlignmentRegion region;
        region.qStart = qStart;
        region.qSpan = qSpan;
        region.tStart = tStart;
        region.tSpan = tSpan;
        region.queryRev = isRev;
        region.type = RegionType::FRONT;
        region.regionId = numRegions;
        ++numRegions;
        ret.regions.emplace_back(std::move(region));
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
        region.qStart = h1.queryPos;
        region.qSpan = h2.queryPos - h1.queryPos;
        region.tStart = h1.targetPos;
        region.tSpan = h2.targetPos - h1.targetPos;
        region.type = RegionType::GLOBAL;
        region.queryRev = isRev;
        region.regionId = numRegions;

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
        ++numRegions;

        // Add the new region.
        ret.regions.emplace_back(std::move(region));
    }

    // Back chunk.
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
        region.qStart = ret.actualQueryEnd;
        region.qSpan = qExtLen;
        region.tStart = ret.actualTargetEnd;
        region.tSpan = tExtLen;
        region.type = RegionType::BACK;
        region.queryRev = isRev;
        region.regionId = numRegions;
        ++numRegions;
        ret.regions.emplace_back(std::move(region));
    }

    return ret;
}

AlignmentResult AlignSingleRegion(const char* targetSeq, int32_t targetLen, const char* querySeqFwd,
                                  const char* querySeqRev, int32_t queryLen,
                                  AlignerBasePtr& alignerGlobal, AlignerBasePtr& alignerExt,
                                  const AlignmentRegion& region)
{
    const char* querySeqInStrand = region.queryRev ? querySeqRev : querySeqFwd;
    const char* targetSeqInStrand = targetSeq;
    int32_t qStart = region.qStart;
    int32_t tStart = region.tStart;
    const int32_t qSpan = region.qSpan;
    const int32_t tSpan = region.tSpan;

    if (qStart >= queryLen || (qStart + qSpan) > queryLen || tStart >= targetLen ||
        (tStart + tSpan) > targetLen) {
        std::ostringstream oss;
        oss << "AlignmentRegion coordinates out of bounds in AlignRegionsGeneric!";
        oss << " qStart = " << qStart << ", qSpan = " << qSpan << ", queryLen = " << queryLen
            << ", tStart = " << tStart << ", tSpan = " << tSpan << ", targetLen = " << targetLen
            << ", regionId = " << region.regionId;
        throw std::runtime_error(oss.str());
    }
    if (targetSeq == NULL || querySeqFwd == NULL || querySeqRev == NULL) {
        std::ostringstream oss;
        oss << "NULL sequence passed to AlignmentRegion!";
        oss << " qStart = " << qStart << ", qSpan = " << qSpan << ", queryLen = " << queryLen
            << ", tStart = " << tStart << ", tSpan = " << tSpan << ", targetLen = " << targetLen
            << ", regionId = " << region.regionId;
        throw std::runtime_error(oss.str());
    }

    // Prepare the reversed front sequence if required.
    std::string qSubSeq;
    std::string tSubSeq;
    if (region.type == RegionType::FRONT) {
        qSubSeq = std::string(querySeqInStrand + region.qStart, region.qSpan);
        tSubSeq = std::string(targetSeqInStrand + region.tStart, region.tSpan);
        std::reverse(qSubSeq.begin(), qSubSeq.end());
        std::reverse(tSubSeq.begin(), tSubSeq.end());
        querySeqInStrand = qSubSeq.c_str();
        targetSeqInStrand = tSubSeq.c_str();
        qStart = 0;
        tStart = 0;
    }

    // Align.
    AlignmentResult alnRes;
    if (region.type == RegionType::FRONT || region.type == RegionType::BACK) {
        alnRes =
            alignerExt->Extend(querySeqInStrand + qStart, qSpan, targetSeqInStrand + tStart, tSpan);
    } else {
        alnRes = alignerGlobal->Global(querySeqInStrand + qStart, qSpan, targetSeqInStrand + tStart,
                                       tSpan);
    }

    if (region.type == RegionType::FRONT) {
        std::reverse(alnRes.cigar.begin(), alnRes.cigar.end());
    }

    return alnRes;
}

RegionsToAlignResults AlignRegionsGeneric(const RegionsToAlign& regions,
                                          AlignerBasePtr& alignerGlobal, AlignerBasePtr& alignerExt)
{
    RegionsToAlignResults ret;

    ret.offsetFrontQuery = 0;
    ret.offsetFrontTarget = 0;
    ret.offsetBackQuery = 0;
    ret.offsetBackTarget = 0;

    for (size_t i = 0; i < regions.regions.size(); ++i) {
        const auto& region = regions.regions[i];

        auto alnRes = AlignSingleRegion(regions.targetSeq, regions.targetLen, regions.querySeqFwd,
                                        regions.querySeqRev, regions.queryLen, alignerGlobal,
                                        alignerExt, region);

        if (region.type == RegionType::FRONT) {
            ret.offsetFrontQuery = alnRes.lastQueryPos;
            ret.offsetFrontTarget = alnRes.lastTargetPos;
        } else if (region.type == RegionType::BACK) {
            ret.offsetBackQuery = alnRes.lastQueryPos;
            ret.offsetBackTarget = alnRes.lastTargetPos;
        }

        // Store the results.
        AlignedRegion alignedRegion{std::move(alnRes.cigar), alnRes.lastQueryPos,
                                    alnRes.lastTargetPos, region.regionId};
        ret.alignedRegions.emplace_back(std::move(alignedRegion));

#ifdef DEBUG_ALIGNMENT_SEEDED
        std::cerr << "[aln region i = " << i << " / " << regions.regions.size() << "]"
                  << " region.regionId = " << region.regionId
                  << ", region.qStart = " << region.qStart
                  << ", region.qEnd = " << (region.qStart + region.qSpan)
                  << ", region.qSpan = " << region.qSpan << ", region.tStart = " << region.tStart
                  << ", region.tEnd = " << (region.tStart + region.tSpan)
                  << ", region.tSpan = " << region.tSpan
                  << ", region.type = " << RegionTypeToString(region.type)
                  << ", region.queryRev = " << (region.queryRev ? "true" : "false")
                  << ", qStart = " << qStart << ", qSpan = " << qSpan << ", tStart = " << tStart
                  << ", tSpan = " << tSpan << "\n"
                  << ", CIGAR: " << ret.alignedRegions.back().cigar.ToStdString() << "\n"
                  << alnRes << "\n\n";

// ValidateCigar(querySeqInStrand + qStart, qSpan, targetSeqInStrand + region.tStart,
//             tSpan, ret.alignedRegions.back().cigar, "Chunk validation.");
#endif
    }

    return ret;
}

OverlapPtr AlignmentSeeded(const OverlapPtr& ovl, const std::vector<SeedHit>& sortedHits,
                           const char* targetSeq, const int32_t targetLen, const char* queryFwd,
                           const char* queryRev, const int32_t queryLen, int32_t minAlignmentSpan,
                           int32_t maxFlankExtensionDist, AlignerBasePtr& alignerGlobal,
                           AlignerBasePtr& alignerExt)
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
    for (const auto& alnRegion : alns.alignedRegions) {
        const auto& currCigar = alnRegion.cigar;
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
