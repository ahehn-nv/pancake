// Author: Ivan Sovic

#ifndef PANCAKE_ALIGNER_BATCH_GPU_H
#define PANCAKE_ALIGNER_BATCH_GPU_H

#include <pacbio/pancake/AlignerBase.h>
#include <pacbio/pancake/AlignerBatchCPU.h>
#include <pacbio/pancake/AlignerFactory.h>
#include <pacbio/pancake/Range.h>
#include <pbbam/Cigar.h>
#include <array>
#include <claraparabricks/genomeworks/cudaaligner/aligner.hpp>
#include <claraparabricks/genomeworks/cudaaligner/alignment.hpp>
#include <claraparabricks/genomeworks/cudaaligner/cudaaligner.hpp>
#include <cstdint>
#include <memory>
#include <vector>

namespace PacBio {
namespace Pancake {

// enum class StatusAddSequencePair
// {
//     OK,
//     SEQUENCE_LEN_BELOW_ZERO,
//     EXCEEDED_MAX_ALIGNMENTS,
//     SEQUENCE_TOO_LONG,
// };

// external enum class StatusAddSequencePair;

int64_t ComputeMaxGPUMemory(int64_t cudaalignerBatches, double maxGPUMemoryFraction);

class AlignerBatchGPU
{
public:
    AlignerBatchGPU(const AlignmentParameters& alnParams, uint32_t maxBandwidth, uint32_t deviceId,
                    double maxGPUMemoryFraction, int64_t maxGPUMemoryCap);
    ~AlignerBatchGPU();

    void Clear();

    StatusAddSequencePair AddSequencePair(const char* query, int32_t queryLen, const char* target,
                                          int32_t targetLen);

    void AlignAll();

    const std::vector<AlignmentResult>& GetAlnResults() const { return alnResults_; }

private:
    AlignmentParameters alnParams_;
    std::unique_ptr<claraparabricks::genomeworks::cudaaligner::Aligner> aligner_;
    cudaStream_t stream_;
    int64_t cudaalignerBatches_;
    int64_t numAlignments_;
    std::vector<int32_t> querySpans_;
    std::vector<int32_t> targetSpans_;
    std::vector<AlignmentResult> alnResults_;
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_ALIGNER_BATCH_GPU_H
