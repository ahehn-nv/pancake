// Authors: Ivan Sovic

#include <pacbio/alignment/AlignmentTools.h>
#include <pacbio/pancake/AlignerBatchGPU.h>
#include <pacbio/util/Util.h>
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

AlignerBatchGPU::AlignerBatchGPU(int32_t numThreads, const AlignmentParameters& alnParams,
                                 uint32_t maxBandwidth, uint32_t deviceId, int64_t maxGPUMemoryCap)
    : AlignerBatchGPU(nullptr, alnParams, maxBandwidth, deviceId, maxGPUMemoryCap)
{
    fafFallback_ = std::make_unique<Parallel::FireAndForget>(numThreads);
    faf_ = fafFallback_.get();
}

AlignerBatchGPU::AlignerBatchGPU(Parallel::FireAndForget* faf, const AlignmentParameters& alnParams,
                                 uint32_t maxBandwidth, uint32_t deviceId, int64_t maxGPUMemoryCap)
    : faf_{faf}, fafFallback_{nullptr}, alnParams_(alnParams)
{
    gpuStream_ = claraparabricks::genomeworks::make_cuda_stream();
    aligner_ = claraparabricks::genomeworks::cudaaligner::create_aligner(
        claraparabricks::genomeworks::cudaaligner::AlignmentType::global_alignment, maxBandwidth,
        gpuStream_.get(), deviceId, maxGPUMemoryCap);
}

AlignerBatchGPU::~AlignerBatchGPU()
{
    aligner_->reset();
    aligner_.reset();
    if (fafFallback_) {
        fafFallback_->Finalize();
    }
}

void AlignerBatchGPU::Clear()
{
    aligner_->reset();
    alnResults_.clear();
    querySpans_.clear();
    targetSpans_.clear();
}

void AlignerBatchGPU::ResetMaxBandwidth(int32_t maxBandwidth)
{
    aligner_->reset_max_bandwidth(maxBandwidth);
}

StatusAddSequencePair AlignerBatchGPU::AddSequencePair(const char* query, int32_t queryLen,
                                                       const char* target, int32_t targetLen)
{
    if (queryLen == 0 || targetLen == 0) {
        // Add a dummy sequence pair to keep the place for the alignment.
        StatusAddSequencePair rv = AddSequencePair_("A", 1, "A", 1);
        querySpans_.back() = queryLen;
        targetSpans_.back() = targetLen;
        return rv;
    }

    return AddSequencePair_(query, queryLen, target, targetLen);
}

StatusAddSequencePair AlignerBatchGPU::AddSequencePair_(const char* query, int32_t queryLen,
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

    // Convert the results to a form we can use downstream.
    alnResults_.resize(querySpans_.size());

    // Determine how many records should land in each thread, spread roughly evenly.
    const int32_t numRecords = alignments.size();
    const std::vector<std::pair<int32_t, int32_t>> jobsPerThread =
        PacBio::Pancake::DistributeJobLoad<int32_t>(faf_ ? faf_->NumThreads() : 1, numRecords);

    // Run the mapping in parallel.
    const auto Submit = [&jobsPerThread, &alignments, this](int32_t i) {
        const int32_t jobStart = jobsPerThread[i].first;
        const int32_t jobEnd = jobsPerThread[i].second;
        WorkerConstructAlignmentResult(jobStart, jobEnd, querySpans_, targetSpans_, alnParams_,
                                       alignments, alnResults_);
    };
    Parallel::Dispatch(faf_, jobsPerThread.size(), Submit);
}

void AlignerBatchGPU::WorkerConstructAlignmentResult(
    int32_t jobStart, int32_t jobEnd, const std::vector<int32_t>& querySpans,
    const std::vector<int32_t>& targetSpans, const AlignmentParameters& alnParams,
    const std::vector<std::shared_ptr<claraparabricks::genomeworks::cudaaligner::Alignment>>&
        alignments,
    std::vector<AlignmentResult>& alnResults)
{
    for (int32_t jobId = jobStart; jobId < jobEnd; ++jobId) {
        auto& alnRes = alnResults[jobId];

        // Handle the edge cases which Cudaaligner does not handle properly.
        if (querySpans[jobId] == 0 || targetSpans[jobId] == 0) {
            alnRes = EdgeCaseAlignmentResult(querySpans[jobId], targetSpans[jobId],
                                             alnParams.matchScore, alnParams.mismatchPenalty,
                                             alnParams.gapOpen1, alnParams.gapExtend1);
            continue;
        }

        alnRes.cigar = CudaalignToCigar_(alignments[jobId]->get_alignment(), alnRes.diffs);
        const int32_t querySpan = alnRes.diffs.numEq + alnRes.diffs.numX + alnRes.diffs.numI;
        const int32_t targetSpan = alnRes.diffs.numEq + alnRes.diffs.numX + alnRes.diffs.numD;
        alnRes.score =
            ScoreCigarAlignment(alnRes.cigar, alnParams.matchScore, alnParams.mismatchPenalty,
                                alnParams.gapOpen1, alnParams.gapExtend1);
        alnRes.valid = (querySpan == querySpans[jobId] && targetSpan == targetSpans[jobId] &&
                        alignments[jobId]->is_optimal());
        alnRes.maxScore = alnRes.score;
        alnRes.zdropped = false;
        alnRes.lastQueryPos = querySpans[jobId];
        alnRes.lastTargetPos = targetSpans[jobId];
        alnRes.maxQueryPos = querySpans[jobId];
        alnRes.maxTargetPos = targetSpans[jobId];
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
