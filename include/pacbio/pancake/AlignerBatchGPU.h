// Author: Ivan Sovic

#ifndef PANCAKE_ALIGNER_BATCH_GPU_H
#define PANCAKE_ALIGNER_BATCH_GPU_H

#include <pacbio/pancake/AlignerBatchCPU.h>
#include <pacbio/pancake/AlignerFactory.h>
#include <pacbio/pancake/Range.h>
#include <pbbam/Cigar.h>
#include <claraparabricks/genomeworks/cudaaligner/aligner.hpp>
#include <claraparabricks/genomeworks/cudaaligner/alignment.hpp>
#include <claraparabricks/genomeworks/cudaaligner/cudaaligner.hpp>
#include <claraparabricks/genomeworks/utils/cudautils.hpp>
#include <cstdint>
#include <memory>
#include <vector>

namespace PacBio {
namespace Pancake {

int64_t ComputeMaxGPUMemory(int64_t cudaalignerBatches, double maxGPUMemoryFraction);
int64_t GetMaxDeviceMemory(int64_t maxDeviceMemoryAllocatorCachingSize);

class AlignerBatchGPU
{
public:
    AlignerBatchGPU(const AlignmentParameters& alnParams, int32_t maxSeqLen, int32_t maxAlignments,
                    uint32_t deviceId, int64_t maxGPUMemoryCap);

    ~AlignerBatchGPU();

    /*
     * Clears the internal states (sequences for alignment and results).
    */
    void Clear();

    /*
     * Reset the cuda aligner to tthe provided maxiumum bandwidth.
    */
    void ResetMaxBandwidth(int32_t maxBandwidth);

    /*
     * Adds a single sequence pair for alignment to the internal state (modifies the state),
     * but does not align the sequences right away. Alignment will be performed only when the
     * AlignAll function is called.
     * The data will be copied internally, so the query and target pointers do not have to
     * remain valid after this function has been called.
     * The return value is either StatusAddSequencePair::OK if the sequences were added successfully,
     * or another value describing the reason for rejecting the sequence pair.
     * Calling this function will clear the internal alignment results.
    */
    StatusAddSequencePair AddSequencePair(const char* query, int32_t queryLen, const char* target,
                                          int32_t targetLen);

    /*
     * Aligns all the sequence pairs added to the aligner, in parallel.
     * This function modifies the internal state, because the alignment results are stored internally.
    */
    void AlignAll();

    /*
     * Const getter for the alignment results. Returns an empty vector if the sequences have not
     * been aligned yet.
    */
    const std::vector<AlignmentResult>& GetAlnResults() const { return alnResults_; }

    /*
     * Non-const getter for the alignment results. Returns an empty vector if the sequences have not
     * been aligned yet.
    */
    std::vector<AlignmentResult>& GetAlnResults() { return alnResults_; }

    size_t BatchSize() const { return querySpans_.size(); }

private:
    AlignmentParameters alnParams_;
    std::unique_ptr<claraparabricks::genomeworks::cudaaligner::Aligner> aligner_;
    claraparabricks::genomeworks::CudaStream gpuStream_;
    std::vector<int32_t> querySpans_;
    std::vector<int32_t> targetSpans_;
    std::vector<AlignmentResult> alnResults_;
    claraparabricks::genomeworks::DefaultDeviceAllocator allocator_;

    StatusAddSequencePair AddSequencePair_(const char* query, int32_t queryLen, const char* target,
                                           int32_t targetLen);

    /*
     * Decodes a single alignment state from the Cudaaligner format to PacBio::Data::CigarOperationType.
    */
    static PacBio::Data::CigarOperationType CudaalignStateToPbbamState_(
        const claraparabricks::genomeworks::cudaaligner::AlignmentState& s);

    /*
     * Helper conversion function to convert the Cudaaligner's alignment into the
     * PacBio::Data::Cigar format. IT also calcualtes the total number of diffs on the fly.
    */
    static PacBio::Data::Cigar CudaalignToCigar_(
        const std::vector<claraparabricks::genomeworks::cudaaligner::AlignmentState>& alignment,
        Alignment::DiffCounts& retDiffs);
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_ALIGNER_BATCH_GPU_H
