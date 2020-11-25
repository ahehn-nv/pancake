// Author: Ivan Sovic

#ifndef PANCAKE_ALIGNER_BATCH_CPU_H
#define PANCAKE_ALIGNER_BATCH_CPU_H

#include <pacbio/pancake/AlignerBase.h>
#include <pacbio/pancake/AlignerFactory.h>
#include <pacbio/pancake/Range.h>
#include <pbbam/Cigar.h>
#include <array>
#include <cstdint>
#include <memory>
#include <vector>

namespace PacBio {
namespace Pancake {

enum class StatusAddSequencePair
{
    OK,
    SEQUENCE_LEN_BELOW_ZERO,
    EXCEEDED_MAX_ALIGNMENTS,
    SEQUENCE_TOO_LONG,
};

class AlignerBatchCPU
{
public:
    AlignerBatchCPU(const AlignerType alnTypeGlobal, const AlignmentParameters& alnParamsGlobal,
                    const AlignerType alnTypeExt, const AlignmentParameters& alnParamsExt);
    ~AlignerBatchCPU();

    void Clear();

    StatusAddSequencePair AddSequencePair(const char* query, int32_t queryLen, const char* target,
                                          int32_t targetLen, bool isGlobalAlignment);
    void AlignAll(int32_t numThreads);

    const std::vector<AlignmentResult>& GetAlnResults() const { return alnResults_; }

private:
    AlignerType alnTypeGlobal_;
    AlignmentParameters alnParamsGlobal_;
    AlignerType alnTypeExt_;
    AlignmentParameters alnParamsExt_;

    // Stores the query and target sequences.
    std::vector<char> bases_;
    int64_t numBases_;

    // Stores locations of query and target sequences.
    std::vector<Range> seqRangesQuery_;
    std::vector<Range> seqRangesTarget_;
    std::vector<bool> isGlobalAlignment_;
    std::vector<AlignmentResult> alnResults_;

    static void Worker_(const std::vector<char>& bases, const std::vector<Range>& seqRangesQuery,
                        const std::vector<Range>& seqRangesTarget,
                        const std::vector<bool>& isGlobalAlignment, int32_t alnStartId,
                        int32_t alnEndId, AlignerBasePtr& alignerGlobal, AlignerBasePtr& alignerExt,
                        std::vector<AlignmentResult>& alnResults);
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_ALIGNER_BATCH_CPU_H
