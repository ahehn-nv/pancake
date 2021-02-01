// Authors: Ivan Sovic

#include <pacbio/alignment/AlignmentTools.h>
#include <pacbio/pancake/AlignerBatchGPU.h>
#include <array>
#include <cstdint>
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
    : alnParams_(alnParams), aligner_(nullptr), stream_(0), cudaalignerBatches_(1)
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
    aligner_.reset();
    alnResults_.clear();
    querySpans_.clear();
    targetSpans_.clear();
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
        return StatusAddSequencePair::EXCEEDED_MAX_ALIGNMENTS;

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
        querySpans_.emplace_back(queryLen);
        targetSpans_.emplace_back(targetLen);
        alnResults_.clear();
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

    // Number of alignments should be the same as number of sequences added.
    if (querySpans_.size() != alignments.size()) {
        throw std::runtime_error(
            "Number of alignments doesn't match number of input sequences, in " +
            std::string(__func__) + ".");
    }

    alnResults_.resize(querySpans_.size());

    for (size_t i = 0; i < alignments.size(); i++) {
        auto& alnRes = alnResults_[i];

        Alignment::DiffCounts diffs;
        alnRes.cigar = CudaalignToCigar_(alignments[i]->get_alignment(), diffs);

        // const int32_t querySpan = diffs.numEq + diffs.numX + diffs.numI;
        // const int32_t targetSpan = diffs.numEq + diffs.numX + diffs.numD;
        alnRes.score =
            ScoreCigarAlignment(alnRes.cigar, alnParams_.matchScore, alnParams_.mismatchPenalty,
                                alnParams_.gapOpen1, alnParams_.gapExtend1);
        // alnRes.valid = (querySpan == querySpans_[i] && targetSpan == targetSpans_[i] && alignments[i]->is_optimal());
        alnRes.valid = (alignments[i]->is_optimal());
        alnRes.maxScore = alnRes.score;
        alnRes.zdropped = false;
        alnRes.lastQueryPos = querySpans_[i];
        alnRes.lastTargetPos = targetSpans_[i];
        alnRes.maxQueryPos = querySpans_[i];
        alnRes.maxTargetPos = targetSpans_[i];
    }
}

PacBio::Data::CigarOperationType AlignerBatchGPU::CudaalignStateToPbbamState_(
    const claraparabricks::genomeworks::cudaaligner::AlignmentState& s)
{
    switch (s) {
        case claraparabricks::genomeworks::cudaaligner::AlignmentState::match:
            return PacBio::Data::CigarOperationType::SEQUENCE_MATCH;
        case claraparabricks::genomeworks::cudaaligner::AlignmentState::mismatch:
            return PacBio::Data::CigarOperationType::SEQUENCE_MISMATCH;
        case claraparabricks::genomeworks::cudaaligner::AlignmentState::insertion:
            return PacBio::Data::CigarOperationType::INSERTION;
        case claraparabricks::genomeworks::cudaaligner::AlignmentState::deletion:
            return PacBio::Data::CigarOperationType::DELETION;
        default:
            return PacBio::Data::CigarOperationType::UNKNOWN_OP;
    }
    return PacBio::Data::CigarOperationType::UNKNOWN_OP;
}

PacBio::Data::Cigar AlignerBatchGPU::CudaalignToCigar_(
    const std::vector<claraparabricks::genomeworks::cudaaligner::AlignmentState>& alignment,
    Alignment::DiffCounts& retDiffs)
{
    retDiffs.Clear();

    if (alignment.empty()) {
        return {};
    }

    std::array<uint32_t, 4> counts{0, 0, 0, 0};

    PacBio::Data::Cigar cigar;
    auto lastState = alignment[0];
    uint32_t count = 0;
    for (auto const& currState : alignment) {
        if (currState == lastState) {
            ++count;
        } else {
            PacBio::Data::CigarOperationType cigarOpType = CudaalignStateToPbbamState_(lastState);
            PacBio::Data::CigarOperation newOp(cigarOpType, count);
            cigar.emplace_back(std::move(newOp));
            counts[lastState] += count;
            count = 1;
            lastState = currState;
        }
    }
    PacBio::Data::CigarOperationType cigarOpType = CudaalignStateToPbbamState_(lastState);
    PacBio::Data::CigarOperation newOp(cigarOpType, count);
    cigar.emplace_back(std::move(newOp));
    counts[lastState] += count;
    retDiffs.numEq = counts[claraparabricks::genomeworks::cudaaligner::AlignmentState::match];
    retDiffs.numX = counts[claraparabricks::genomeworks::cudaaligner::AlignmentState::mismatch];
    retDiffs.numI = counts[claraparabricks::genomeworks::cudaaligner::AlignmentState::insertion];
    retDiffs.numD = counts[claraparabricks::genomeworks::cudaaligner::AlignmentState::deletion];

    return cigar;
}

}  // namespace Pancake
}  // namespace PacBio
