// Authors: Ivan Sovic

#include <pacbio/alignment/AlignmentTools.h>
#include <pacbio/pancake/AlignerBatchGPU.h>
#include <cstring>
#include <iostream>
#include <sstream>

namespace PacBio {
namespace Pancake {

int64_t ComputeMaxGPUMemory(int64_t cudaalignerBatches, double maxGPUMemoryFraction)
{
    size_t freeMem = 0, totalMem = 0;

    GW_CU_CHECK_ERR(cudaMemGetInfo(&freeMem, &totalMem));
    const size_t freeUsableMemory = static_cast<double>(freeMem) * maxGPUMemoryFraction;
    const int64_t usableMemoryPerAligner = freeUsableMemory / cudaalignerBatches;

    std::cerr << "cudaalignerBatches = " << cudaalignerBatches << "\n";
    std::cerr << "freeMem = " << freeMem / (1024.0 * 1024.0) << " MB\n";
    std::cerr << "totalMem = " << totalMem / (1024.0 * 1024.0) << " MB\n";
    std::cerr << "usableMemoryPerAligner = " << usableMemoryPerAligner / (1024.0 * 1024.0)
              << " MB \n";

    return usableMemoryPerAligner;
}

AlignerBatchGPU::AlignerBatchGPU(const AlignmentParameters& alnParams, uint32_t maxBandwidth,
                                 uint32_t deviceId, double maxGPUMemoryFraction,
                                 int64_t maxGPUMemoryCap)
    : alnParams_(alnParams)
    , aligner_(nullptr)
    , stream_(0)
    , cudaalignerBatches_(1)
    , numAlignments_(0)
{
    GW_CU_CHECK_ERR(cudaSetDevice(deviceId));
    GW_CU_CHECK_ERR(cudaStreamCreate(&stream_));

    const int64_t requestedMemory =
        std::min(maxGPUMemoryCap, PacBio::Pancake::ComputeMaxGPUMemory(1, maxGPUMemoryFraction));

    aligner_ = claraparabricks::genomeworks::cudaaligner::create_aligner(
        claraparabricks::genomeworks::cudaaligner::AlignmentType::global_alignment, maxBandwidth,
        stream_, deviceId, requestedMemory);
}

AlignerBatchGPU::~AlignerBatchGPU()
{
    aligner_.reset();
    GW_CU_CHECK_ERR(cudaStreamDestroy(stream_));
}

void AlignerBatchGPU::Clear()
{
    alnResults_.clear();
    alnResults_.resize(numAlignments_);
    numAlignments_ = 0;
}

StatusAddSequencePair AlignerBatchGPU::AddSequencePair(const char* query, int32_t queryLen,
                                                       const char* target, int32_t targetLen)
{
    if (queryLen < 0 || targetLen < 0) {
        return StatusAddSequencePair::SEQUENCE_LEN_BELOW_ZERO;
    }

    const claraparabricks::genomeworks::cudaaligner::StatusType s =
        aligner_->add_alignment(target, targetLen, query, queryLen);

    if (s == claraparabricks::genomeworks::cudaaligner::StatusType::exceeded_max_alignments) {
        throw std::runtime_error("Exceeded max alignments.");

    } else if (s == claraparabricks::genomeworks::cudaaligner::StatusType::
                        exceeded_max_alignment_difference ||
               s == claraparabricks::genomeworks::cudaaligner::StatusType::exceeded_max_length) {
        // Do nothing as this case will be handled by CPU aligner.
        throw std::runtime_error(
            "Not implemented yet: s == StatusType::exceeded_max_alignment_difference || s == "
            "StatusType::exceeded_max_length.\n");

    } else if (s != claraparabricks::genomeworks::cudaaligner::StatusType::success) {
        // fprintf(stderr, "Unknown error in cuda aligner!\n");
        throw std::runtime_error("Unknown error in cuda aligner!\n");

    } else {
        ++numAlignments_;
        querySpans_.emplace_back(queryLen);
        targetSpans_.emplace_back(targetLen);
        // std::cerr << "Added:\n  -> Q: '" << std::string(query, queryLen) << "'\n  -> T: '" << std::string(target, targetLen) << "'\n";
    }

    return StatusAddSequencePair::OK;
}

void AlignerBatchGPU::AlignAll()
{
    alnResults_.clear();

    aligner_->align_all();
    aligner_->sync_alignments();

    const std::vector<std::shared_ptr<claraparabricks::genomeworks::cudaaligner::Alignment>>&
        alignments = aligner_->get_alignments();

    // Number of alignments should be the same as number of overlaps.
    if (numAlignments_ != static_cast<int64_t>(alignments.size())) {
        throw std::runtime_error(
            "Number of alignments doesn't match number of overlaps in cudaaligner.");
    }

    alnResults_.resize(numAlignments_);

    for (size_t i = 0; i < alignments.size(); i++) {
        auto& alnRes = alnResults_[i];

        const std::string cigar = alignments[i]->convert_to_cigar(
            claraparabricks::genomeworks::cudaaligner::CigarFormat::extended);

        // std::cerr << "[i = " << i << "] CUDA CIGAR: " << cigar << "\n";

        alnRes.cigar = PacBio::BAM::Cigar(cigar);

        const Alignment::DiffCounts diffs = CigarDiffCounts(cigar);
        const int32_t querySpan = diffs.numEq + diffs.numX + diffs.numI;
        const int32_t targetSpan = diffs.numEq + diffs.numX + diffs.numD;
        alnRes.score =
            ScoreCigarAlignment(alnRes.cigar, alnParams_.matchScore, alnParams_.mismatchPenalty,
                                alnParams_.gapOpen1, alnParams_.gapExtend1);
        alnRes.valid = (querySpan == querySpans_[i] && targetSpan == targetSpans_[i]);
        alnRes.maxScore = alnRes.score;
        alnRes.zdropped = false;
        alnRes.lastQueryPos = querySpans_[i];
        alnRes.lastTargetPos = targetSpans_[i];
        alnRes.maxQueryPos = querySpans_[i];
        alnRes.maxTargetPos = targetSpans_[i];
    }
}

}  // namespace Pancake
}  // namespace PacBio
