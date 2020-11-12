// Authors: Ivan Sovic

#include <pacbio/pancake/AlignerKSW2.h>
#include <pacbio/pancake/Lookups.h>
#include <cstring>

namespace PacBio {
namespace Pancake {

mm_tbuf_t* mm_tbuf_init(void)
{
    mm_tbuf_t* b;
    b = (mm_tbuf_t*)calloc(1, sizeof(mm_tbuf_t));
    // if (!(mm_dbg_flag & 1)) b->km = km_init();              // This flag means: #define MM_DBG_NO_KALLOC     0x1
    return b;
}

void mm_tbuf_destroy(mm_tbuf_t* b)
{
    if (b == 0) return;
    km_destroy(b->km);
    free(b);
}

std::shared_ptr<AlignerBase> CreateAlignerKSW2(const AlignmentParameters& opt)
{
    return std::shared_ptr<AlignerBase>(new AlignerKSW2(opt));
}

AlignerKSW2::AlignerKSW2(const AlignmentParameters& opt)
    : opt_(opt), buffer_(mm_tbuf_init(), mm_tbuf_destroy)
{
    // Minimap2 alignment matrix.
    // From Minimap2: int sc_ambi; // score when one or both bases are "N"
    const int32_t scAmbi = 1;
    ksw_gen_simple_mat_(5, mat_, opt.matchScore, opt.mismatchPenalty, scAmbi);
}

AlignerKSW2::~AlignerKSW2() = default;

AlignmentResult AlignerKSW2::Global(const char* qseq, int64_t qlen, const char* tseq, int64_t tlen)
{
    const int32_t extra_flag = 0;

    // Bandwidth heuristic, as it is in Minimap2.
    const int32_t bw = (int)(opt_.alignBandwidth * 1.5 + 1.);

    // Memory allocations required for KSW2.
    ksw_extz_t ez;
    memset(&ez, 0, sizeof(ksw_extz_t));

    // Convert the subsequence's alphabet from ACTG to [0123].
    const std::vector<uint8_t> qseqInt = ConvertSeqAlphabet_(qseq, qlen, &BaseToTwobit[0]);
    const std::vector<uint8_t> tseqInt = ConvertSeqAlphabet_(tseq, tlen, &BaseToTwobit[0]);

    // Compute the actual bandwidth. If this was a long join, we need to allow more room.
    const int32_t longestSpan = std::max(qlen, tlen);
    const int32_t spanDiff = std::abs(tlen - qlen);

    // This is slighlty modified compared to how it is in Minimap2:
    // const int32_t actualBandwidth =
    //     (h2.CheckFlagLongJoin() || spanDiff > (bw / 2)) ? longestSpan : bw;
    const int32_t actualBandwidth = (spanDiff > (bw / 2)) ? longestSpan : bw;

    // First pass: with approximate Z-drop
    mm_align_pair_(buffer_->km, qlen, &qseqInt[0], tlen, &tseqInt[0], mat_, actualBandwidth, -1, -1,
                   extra_flag | KSW_EZ_APPROX_MAX, &ez, opt_.gapOpen1, opt_.gapExtend1,
                   opt_.gapOpen2, opt_.gapExtend2);

    // // First pass: with approximate Z-drop
    // mm_align_pair_(buffer_->km, &mopt_, qlen, &qseqInt[0], tlen, &tseqInt[0], mat_, actualBandwidth,
    //               -1, mopt_.zdrop, extra_flag | KSW_EZ_APPROX_MAX, &ez);

    // // Test Z-drop and inversion Z-drop
    // const int32_t zdrop_code =
    //     mm_test_zdrop_(buffer_->km, &mopt_, &qseqInt[0], &tseqInt[0], ez.n_cigar, ez.cigar, mat_);
    // if (zdrop_code != 0) {
    //     // second pass: lift approximate
    //     mm_align_pair_(buffer_->km, &mopt_, qlen, &qseqInt[0], tlen, &tseqInt[0], mat_,
    //                   actualBandwidth, -1, zdrop_code == 2 ? mopt_.zdrop_inv : mopt_.zdrop,
    //                   extra_flag, &ez);
    // }
    auto currCigar = ConvertMinimap2CigarToPbbam_(ez.cigar, ez.n_cigar, qseqInt, tseqInt);

    AlignmentResult ret;
    ret.cigar = std::move(currCigar);
    // ret.valid = ez.reach_end;
    ret.valid = true;
    ret.lastQueryPos = qlen;
    ret.lastTargetPos = tlen;
    // ret.lastQueryPos = (ez.reach_end ? qlen : ez.max_q + 1);
    // ret.lastTargetPos = (ez.reach_end ? ez.mqe_t + 1 : ez.max_t + 1);
    ret.maxQueryPos = ez.max_q;
    ret.maxTargetPos = ez.max_t;
    ret.score = ez.score;
    ret.maxScore = ez.max;
    ret.zdropped = ez.zdropped;

    // Free KSW2 memory.
    kfree(buffer_->km, ez.cigar);

    return ret;
}

AlignmentResult AlignerKSW2::Extend(const char* qseq, int64_t qlen, const char* tseq, int64_t tlen)
{
    const int32_t extra_flag = 0;

    // Bandwidth heuristic, as it is in Minimap2.
    const int32_t bw = (int)(opt_.alignBandwidth * 1.5 + 1.);

    // Memory allocations required for KSW2.
    ksw_extz_t ez;
    memset(&ez, 0, sizeof(ksw_extz_t));

    // Convert the subsequence's alphabet from ACTG to [0123].
    const std::vector<uint8_t> qseqInt = ConvertSeqAlphabet_(qseq, qlen, &BaseToTwobit[0]);
    const std::vector<uint8_t> tseqInt = ConvertSeqAlphabet_(tseq, tlen, &BaseToTwobit[0]);

    mm_align_pair_(buffer_->km, qlen, &qseqInt[0], tlen, &tseqInt[0], mat_, bw, opt_.endBonus,
                   opt_.zdrop, extra_flag | KSW_EZ_EXTZ_ONLY | KSW_EZ_RIGHT, &ez, opt_.gapOpen1,
                   opt_.gapExtend1, opt_.gapOpen2, opt_.gapExtend2);

    auto currCigar = ConvertMinimap2CigarToPbbam_(ez.cigar, ez.n_cigar, qseqInt, tseqInt);

    AlignmentResult ret;
    ret.cigar = std::move(currCigar);
    ret.valid = ez.reach_end;
    ret.lastQueryPos = (ez.reach_end ? qlen : ez.max_q + 1);
    ret.lastTargetPos = (ez.reach_end ? ez.mqe_t + 1 : ez.max_t + 1);
    ret.maxQueryPos = ez.max_q;
    ret.maxTargetPos = ez.max_t;
    ret.score = ez.score;
    ret.maxScore = ez.max;
    ret.zdropped = ez.zdropped;

    // Free KSW2 memory.
    kfree(buffer_->km, ez.cigar);

    return ret;
}

std::vector<uint8_t> AlignerKSW2::ConvertSeqAlphabet_(const char* seq, size_t seqlen,
                                                      const int8_t* conv_table)
{
    std::vector<uint8_t> ret(seqlen);
    for (size_t i = 0; i < seqlen; i++) {
        ret[i] = (int8_t)conv_table[(uint8_t)seq[i]];
    }
    return ret;
}

PacBio::Data::Cigar AlignerKSW2::ConvertMinimap2CigarToPbbam_(uint32_t* mm2Cigar, int32_t cigarLen,
                                                              const std::vector<uint8_t>& qseq,
                                                              const std::vector<uint8_t>& tseq)
{
    PacBio::Data::Cigar cigar;

    int32_t qPos = 0;
    int32_t tPos = 0;
    for (int32_t opId = 0; opId < cigarLen; ++opId) {
        char op = "MIDNSH"[mm2Cigar[opId] & 0xf];
        int32_t count = mm2Cigar[opId] >> 4;
        // cigar.emplace_back(Data::CigarOperation(op, count));

        if (op == 'M') {
            // std::cerr << "[opId = " << opId << "] op: " << count << op << "\n";

            char prevOp = 'u';
            int32_t span = 0;
            for (int32_t m = 0; m < count; ++m) {
                char currOp = 'X';
                if (qseq[qPos + m] == tseq[tPos + m]) {
                    currOp = '=';
                }
                if (currOp == prevOp) {
                    ++span;
                } else {
                    if (span > 0) {
                        cigar.emplace_back(Data::CigarOperation(prevOp, span));
                        // std::cerr << "  -> Added (mid):  " << span << prevOp << "\n";
                    }
                    span = 1;
                }
                prevOp = currOp;
                // ++qPos;
                // ++tPos;
            }
            if (span > 0) {
                cigar.emplace_back(Data::CigarOperation(prevOp, span));
                // std::cerr << "  -> Added (last): " << span << prevOp << "\n";
            }
            qPos += count;
            tPos += count;
        } else if (op == '=' || op == 'X') {
            cigar.emplace_back(Data::CigarOperation(op, count));
            qPos += count;
            tPos += count;
        } else if (op == 'I' || op == 'S') {
            cigar.emplace_back(Data::CigarOperation(op, count));
            qPos += count;
        } else if (op == 'D' || op == 'N') {
            cigar.emplace_back(Data::CigarOperation(op, count));
            tPos += count;
        }
    }
    return cigar;
}

void AlignerKSW2::mm_align_pair_(void* km, int qlen, const uint8_t* qseq, int tlen,
                                 const uint8_t* tseq, const int8_t* mat, int w, int end_bonus,
                                 int zdrop, int flag, ksw_extz_t* ez, int q, int e, int q2, int e2)
{
    // if (mm_dbg_flag & MM_DBG_PRINT_ALN_SEQ) {
    //     int i;
    //     fprintf(stderr, "===> q=(%d,%d), e=(%d,%d), bw=%d, flag=%d, zdrop=%d <===\n", opt->q,
    //             opt->q2, opt->e, opt->e2, w, flag, opt->zdrop);
    //     for (i = 0; i < tlen; ++i)
    //         fputc("ACGTN"[tseq[i]], stderr);
    //     fputc('\n', stderr);
    //     for (i = 0; i < qlen; ++i)
    //         fputc("ACGTN"[qseq[i]], stderr);
    //     fputc('\n', stderr);
    // }
    // if (opt->max_sw_mat > 0 && (int64_t)tlen * qlen > opt->max_sw_mat) {
    //     ksw_reset_extz(ez);
    //     ez->zdropped = 1;
    // }
    // else if (opt->flag & MM_F_SPLICE)
    //     ksw_exts2_sse(km, qlen, qseq, tlen, tseq, 5, mat, opt->q, opt->e, opt->q2, opt->noncan,
    //                   zdrop, flag, ez);
    if (q == q2 && e == e2)
        ksw_extz2_sse(km, qlen, qseq, tlen, tseq, 5, mat, q, e, w, zdrop, end_bonus, flag, ez);
    else
        ksw_extd2_sse(km, qlen, qseq, tlen, tseq, 5, mat, q, e, q2, e2, w, zdrop, end_bonus, flag,
                      ez);
    // if (mm_dbg_flag & MM_DBG_PRINT_ALN_SEQ) {
    //     int i;
    //     fprintf(stderr, "score=%d, cigar=", ez->score);
    //     for (i = 0; i < ez->n_cigar; ++i)
    //         fprintf(stderr, "%d%c", ez->cigar[i] >> 4, "MIDN"[ez->cigar[i] & 0xf]);
    //     fprintf(stderr, "\n");
    // }
}

// inline void update_max_zdrop_(int32_t score, int i, int j, int32_t* max, int* max_i, int* max_j,
//                               int e, int* max_zdrop, int pos[2][2])
// {
//     if (score < *max) {
//         int li = i - *max_i;
//         int lj = j - *max_j;
//         int diff = li > lj ? li - lj : lj - li;
//         int z = *max - score - diff * e;
//         if (z > *max_zdrop) {
//             *max_zdrop = z;
//             pos[0][0] = *max_i, pos[0][1] = i + 1;
//             pos[1][0] = *max_j, pos[1][1] = j + 1;
//         }
//     } else
//         *max = score, *max_i = i, *max_j = j;
// }

// int AlignerKSW2::mm_test_zdrop_(void* km, const mm_mapopt_t* opt, const uint8_t* qseq,
//                                 const uint8_t* tseq, uint32_t n_cigar, uint32_t* cigar,
//                                 const int8_t* mat)
// {
//     uint32_t k;
//     int32_t score = 0, max = INT32_MIN, max_i = -1, max_j = -1, i = 0, j = 0, max_zdrop = 0;
//     int pos[2][2] = {{-1, -1}, {-1, -1}}, q_len, t_len;

//     // find the score and the region where score drops most along diagonal
//     for (k = 0, score = 0; k < n_cigar; ++k) {
//         uint32_t l, op = cigar[k] & 0xf, len = cigar[k] >> 4;
//         if (op == 0) {
//             for (l = 0; l < len; ++l) {
//                 score += mat[tseq[i + l] * 5 + qseq[j + l]];
//                 update_max_zdrop_(score, i + l, j + l, &max, &max_i, &max_j, opt->e, &max_zdrop,
//                                   pos);
//             }
//             i += len, j += len;
//         } else if (op == 1 || op == 2 || op == 3) {
//             score -= opt->q + opt->e * len;
//             if (op == 1)
//                 j += len;  // insertion
//             else
//                 i += len;  // deletion
//             update_max_zdrop_(score, i, j, &max, &max_i, &max_j, opt->e, &max_zdrop, pos);
//         }
//     }

//     // test if there is an inversion in the most dropped region
//     q_len = pos[1][1] - pos[1][0], t_len = pos[0][1] - pos[0][0];
//     if (!(opt->flag & (MM_F_SPLICE | MM_F_SR | MM_F_FOR_ONLY | MM_F_REV_ONLY)) &&
//         max_zdrop > opt->zdrop_inv && q_len < opt->max_gap && t_len < opt->max_gap) {
//         uint8_t* qseq2;
//         void* qp;
//         int q_off, t_off;
//         qseq2 = (uint8_t*)kmalloc(km, q_len);
//         for (i = 0; i < q_len; ++i) {
//             int c = qseq[pos[1][1] - i - 1];
//             qseq2[i] = c >= 4 ? 4 : 3 - c;
//         }
//         qp = ksw_ll_qinit(km, 2, q_len, qseq2, 5, mat);
//         score = ksw_ll_i16(qp, t_len, tseq + pos[0][0], opt->q, opt->e, &q_off, &t_off);
//         kfree(km, qseq2);
//         kfree(km, qp);
//         if (score >= opt->min_chain_score * opt->a && score >= opt->min_dp_max)
//             return 2;  // there is a potential inversion
//     }
//     return max_zdrop > opt->zdrop ? 1 : 0;
// }

void AlignerKSW2::ksw_gen_simple_mat_(int m, int8_t* mat, int8_t a, int8_t b, int8_t sc_ambi)
{
    int i, j;
    a = a < 0 ? -a : a;
    b = b > 0 ? -b : b;
    sc_ambi = sc_ambi > 0 ? -sc_ambi : sc_ambi;
    for (i = 0; i < m - 1; ++i) {
        for (j = 0; j < m - 1; ++j)
            mat[i * m + j] = i == j ? a : b;
        mat[i * m + m - 1] = sc_ambi;
    }
    for (j = 0; j < m; ++j)
        mat[(m - 1) * m + j] = sc_ambi;
}

}  // namespace Pancake
}  // namespace PacBio
