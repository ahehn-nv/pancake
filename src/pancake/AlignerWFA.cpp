// Authors: Ivan Sovic

#include <pacbio/alignment/AlignmentTools.h>
#include <pacbio/pancake/AlignerWFA.h>

namespace PacBio {
namespace Pancake {

std::shared_ptr<AlignerBase> CreateAlignerWFA(const AlignmentParameters& opt)
{
    return std::shared_ptr<AlignerBase>(new AlignerWFA(opt));
}

AlignerWFA::AlignerWFA(const AlignmentParameters& opt) : opt_(opt), mmAllocator_(NULL)
{
    mmAllocator_ = mm_allocator_new(BUFFER_SIZE_8M);

    penalties_ = {
        .match = -opt.matchScore,
        .mismatch = std::abs(opt.mismatchPenalty),
        .gap_opening = std::abs(opt.gapOpen1),
        .gap_extension = std::abs(opt.gapExtend1),
    };
}

AlignerWFA::~AlignerWFA() { mm_allocator_delete(mmAllocator_); }

AlignmentResult AlignerWFA::Global(const char* qseq, int64_t qlen, const char* tseq, int64_t tlen)
{
    /*
     * Non-adaptive, exact aligner.
    */
    // // Init Affine-WFA.
    // affine_wavefronts_t* affineWavefronts =
    //     affine_wavefronts_new_complete(qlen, tlen, &penalties_, NULL, mmAllocator_);
    // // Align.
    // affine_wavefronts_align(affineWavefronts, qseq, qlen, tseq, tlen);

    /*
     * Adaptive aligner.
    */
    const int min_wavefront_length = 10;
    const int max_distance_threshold = 50;
    // const int min_wavefront_length = 1000;
    // const int max_distance_threshold = 2000;
    // Init Affine-WFA
    affine_wavefronts_t* affineWavefronts = affine_wavefronts_new_reduced(
        qlen, tlen, &penalties_, min_wavefront_length, max_distance_threshold, NULL, mmAllocator_);
    // Align.
    affine_wavefronts_align(affineWavefronts, qseq, qlen, tseq, tlen);

    auto cigar = WFAAlignmentToCigar_(affineWavefronts->edit_cigar);
    const int32_t score = edit_cigar_score_gap_affine(&affineWavefronts->edit_cigar, &penalties_);
    bool valid = true;

    try {
        cigar = NormalizeCigar(qseq, qlen, tseq, tlen, cigar);
    } catch (std::exception& e) {
        valid = false;
        cigar.clear();
    }

    affine_wavefronts_delete(affineWavefronts);

    // Adjust the results.
    AlignmentResult ret;
    ret.cigar = std::move(cigar);
    ret.score = score;
    ret.valid = valid;
    ret.maxScore = score;
    ret.zdropped = false;
    ret.lastQueryPos = qlen;
    ret.lastTargetPos = tlen;
    ret.maxQueryPos = qlen;
    ret.maxTargetPos = tlen;
    return ret;
}

AlignmentResult AlignerWFA::Extend(const char* qseq, int64_t qlen, const char* tseq, int64_t tlen)
{
    AlignmentResult ret;
    ret.valid = false;
    ret.lastQueryPos = 0;
    ret.lastTargetPos = 0;
    return ret;
}

PacBio::BAM::Cigar AlignerWFA::WFAAlignmentToCigar_(const edit_cigar_t& wfaCigar)
{
    const int32_t alnLen = wfaCigar.end_offset - wfaCigar.begin_offset;
    if (alnLen <= 0) {
        return {};
    }

    // Lookup table for conversion.
    std::array<PacBio::BAM::CigarOperationType, 256> opToCigar;
    for (size_t i = 0; i < opToCigar.size(); ++i) {
        opToCigar[i] = PacBio::BAM::CigarOperationType::UNKNOWN_OP;
    }
    opToCigar['M'] = PacBio::BAM::CigarOperationType::SEQUENCE_MATCH;
    opToCigar['='] = PacBio::BAM::CigarOperationType::SEQUENCE_MATCH;
    opToCigar['X'] = PacBio::BAM::CigarOperationType::SEQUENCE_MISMATCH;
    // Swap insertion and deletion operations because WFA considers these with respect
    // to the "text" (target) sequence.
    opToCigar['I'] = PacBio::BAM::CigarOperationType::DELETION;
    opToCigar['D'] = PacBio::BAM::CigarOperationType::INSERTION;

    PacBio::BAM::CigarOperationType prevOp = PacBio::BAM::CigarOperationType::UNKNOWN_OP;
    int32_t count = 0;
    PacBio::BAM::Cigar ret;
    for (int32_t i = wfaCigar.begin_offset, j = 0; i <= wfaCigar.end_offset; ++i, ++j) {
        const char op = wfaCigar.operations[i];

        if (j == alnLen ||
            (opToCigar[op] != prevOp && prevOp != PacBio::BAM::CigarOperationType::UNKNOWN_OP)) {
            ret.emplace_back(PacBio::BAM::CigarOperation(prevOp, count));
            count = 0;
        }
        if (j < alnLen) {
            prevOp = opToCigar[op];
            count += 1;
        }
    }
    return ret;
}

}  // namespace Pancake
}  // namespace PacBio
