// Author: Ivan Sovic

#ifndef PANCAKE_ALIGNER_WFA_H
#define PANCAKE_ALIGNER_WFA_H

#include <gap_affine/affine_wavefront_align.h>
#include <lib/wfa/system/mm_allocator.h>
#include <pacbio/pancake/AlignerBase.h>
#include <pbbam/Cigar.h>
#include <array>
#include <cstdint>
#include <memory>
#include <vector>

namespace PacBio {
namespace Pancake {

class AlignerWFA;
std::shared_ptr<AlignerBase> CreateAlignerWFA(const AlignmentParameters& opt);

class AlignerWFA : public AlignerBase
{
public:
    AlignerWFA(const AlignmentParameters& opt);
    ~AlignerWFA() override;

    AlignmentResult Global(const char* qseq, int64_t qlen, const char* tseq, int64_t tlen) override;
    AlignmentResult Extend(const char* qseq, int64_t qlen, const char* tseq, int64_t tlen) override;

private:
    AlignmentParameters opt_;
    mm_allocator_t* mmAllocator_;
    affine_penalties_t penalties_;

    static PacBio::BAM::Cigar WFAAlignmentToCigar_(const edit_cigar_t& wfaCigar);
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_ALIGNER_WFA_H
