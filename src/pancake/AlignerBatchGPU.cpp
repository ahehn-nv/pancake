// Copyright (c) 2019, Pacific Biosciences of California, Inc.
// All rights reserved.
// See LICENSE.txt.
//
// Contributions from NVIDIA are Copyright (c) 2021, NVIDIA Corporation.
// All rights reserved.
// SPDX-License-Identifier: BSD-3-Clause-Clear
//
// Authors: Ivan Sovic

#include <pacbio/alignment/AlignmentTools.h>
#include <pacbio/pancake/AlignerBatchGPU.h>
#include <pacbio/pancake/AlignerBatchGPUGenomeWorksInterface.h>
#include <pacbio/util/Util.h>
#include <claraparabricks/genomeworks/utils/device_buffer.hpp>
#include <claraparabricks/genomeworks/utils/pinned_host_vector.hpp>
#include <cstdint>
#include <iostream>
#include <sstream>
#include <limits>

namespace gw = claraparabricks::genomeworks;

namespace PacBio {
namespace Pancake {

struct AlignerBatchGPU::AlignerBatchGPUHostBuffers
{
    claraparabricks::genomeworks::pinned_host_vector<int32_t> starts;
    claraparabricks::genomeworks::pinned_host_vector<int32_t> index;
    claraparabricks::genomeworks::pinned_host_vector<int4> diffs;
    claraparabricks::genomeworks::pinned_host_vector<int64_t> scores;
    claraparabricks::genomeworks::pinned_host_vector<PacBio::Data::CigarOperation> cigarOpsBuffer;

    void clearAndResize(int64_t cigarOpsBufferLen, int32_t numberAlns)
    {
        if (cigarOpsBuffer.size() < static_cast<size_t>(cigarOpsBufferLen)) {
            cigarOpsBuffer.clear();
            cigarOpsBuffer.shrink_to_fit();
        }
        if (starts.size() < static_cast<size_t>(numberAlns)) {
            // add 20% to avoid later reallocations
            const int32_t newSize = numberAlns * 12 / 10;
            starts.clear();
            index.clear();
            diffs.clear();
            scores.clear();
            starts.shrink_to_fit();
            index.shrink_to_fit();
            diffs.shrink_to_fit();
            scores.shrink_to_fit();
            starts.resize(newSize + 1);
            index.resize(newSize);
            diffs.resize(newSize);
            scores.resize(newSize);
        }
        if (cigarOpsBuffer.size() < static_cast<size_t>(cigarOpsBufferLen)) {
            // add 20% to avoid later reallocations
            const int32_t newSize = cigarOpsBufferLen * 12 / 10;
            cigarOpsBuffer.clear();
            cigarOpsBuffer.shrink_to_fit();
            cigarOpsBuffer.resize(newSize);
        }
    }

    void free()
    {
        cigarOpsBuffer.clear();
        cigarOpsBuffer.shrink_to_fit();
        starts.clear();
        index.clear();
        diffs.clear();
        scores.clear();
        starts.shrink_to_fit();
        index.shrink_to_fit();
        diffs.shrink_to_fit();
        scores.shrink_to_fit();
    }
};

static int32_t findLongestCigarLength(const gw::cudaaligner::DeviceAlignmentsPtrs& aln, int32_t numberOfAlignments, AlignerBatchGPU::AlignerBatchGPUHostBuffers* hostBuffers, cudaStream_t stream)
{
    int64_t longestCigarLength = 0;
    gw::cudautils::device_copy_n_async(aln.starts, numberOfAlignments+1, hostBuffers->starts.data(), stream);
    GW_CU_CHECK_ERR(cudaStreamSynchronize(stream));
    for(auto it = begin(hostBuffers->starts) + 1; it < end(hostBuffers->starts); ++it)
    {
        const int64_t cigarLength = *it - *(it - 1);
        assert(cigarLength >= 0);
        longestCigarLength = std::max(longestCigarLength, cigarLength);
    }
    assert(longestCigarLength < std::numeric_limits<int32_t>::max());
    return static_cast<int32_t>(longestCigarLength);
}

void RetrieveResultsAsPacBioCigar(
    AlignerBatchGPU::AlignerBatchGPUHostBuffers* hostBuffers,
    const std::vector<int32_t>& querySpans, const std::vector<int32_t>& targetSpans,
    const claraparabricks::genomeworks::cudaaligner::FixedBandAligner* aligner,
    std::vector<AlignmentResult>& results, int32_t matchScore, int32_t mismatchScore,
    int32_t gapOpenScore, int32_t gapExtScore)
{
    GW_NVTX_RANGE(profiler, "pancake retrieve results");
    namespace gw = claraparabricks::genomeworks;
    if (hostBuffers == nullptr || aligner == nullptr) {
        throw std::runtime_error("hostBuffers or aligner should not be nullptr " +
                                 std::string(__func__) + ".");
    }
    cudaStream_t stream = aligner->get_stream();
    auto allocator = aligner->get_device_allocator();
    GW_CU_CHECK_ERR(cudaStreamSynchronize(stream));
    gw::cudaaligner::DeviceAlignmentsPtrs aln = aligner->get_alignments_device();
    const int32_t numberOfAlignments = aln.n_alignments;

    // Number of alignments should be the same as number of sequences added.
    if (querySpans.size() != static_cast<size_t>(numberOfAlignments)) {
        throw std::runtime_error(
            "Number of alignments doesn't match number of input sequences, in " +
            std::string(__func__) + ".");
    }

    results.clear();

    if (numberOfAlignments == 0) {
        return;
    }

    gw::device_buffer<int4> diffs_device(numberOfAlignments, allocator, stream);
    gw::device_buffer<int64_t> scores_device(numberOfAlignments, allocator, stream);
    gw::device_buffer<PacBio::Data::CigarOperation> pacbio_cigars_device(0, allocator, stream);

    int64_t cigarBufSize       = aln.total_length;
    int32_t longestCigarLength = 0;
    while (1)
    {
        if(cigarBufSize < longestCigarLength)
        {
            throw std::runtime_error("Could not allocate enough device or pinned host memory in " + std::string(__func__) + ".");
        }
        try
        {
            pacbio_cigars_device.clear_and_resize(cigarBufSize);
            hostBuffers->clearAndResize(cigarBufSize, numberOfAlignments);
            break; // allocation successful
        }
        catch (const std::bad_alloc& except)
        {
        }
        catch(const gw::device_memory_allocation_exception& except)
        {
        }

        if (hostBuffers->index.size() < static_cast<size_t>(numberOfAlignments))
        {
            throw std::runtime_error("Could not allocate enough pinned host memory in " + std::string(__func__) + ".");
        }
        hostBuffers->free();

        if (longestCigarLength == 0)
        {
            longestCigarLength = findLongestCigarLength(aln, numberOfAlignments, hostBuffers, stream);
        }
        cigarBufSize /= 2;
    }

    gw::cudautils::device_copy_n_async(aln.starts, numberOfAlignments + 1, hostBuffers->starts.data(), stream);
    gw::cudautils::device_copy_n_async(aln.lengths, numberOfAlignments, hostBuffers->index.data(), stream);
    results.resize(numberOfAlignments);
    GW_CU_CHECK_ERR(cudaStreamSynchronize(stream));
    int32_t completedAlignments = 0;
    while (completedAlignments < numberOfAlignments)
    {
        GWInterface::RunConvertToPacBioCigarAndScoreGpuKernel(pacbio_cigars_device.data(), diffs_device.data(), scores_device.data(), aln, cigarBufSize, completedAlignments, matchScore, mismatchScore, gapOpenScore, gapExtScore, stream, aligner->get_device());

        gw::cudautils::device_copy_n_async(diffs_device.data(), numberOfAlignments, hostBuffers->diffs.data(), stream);
        gw::cudautils::device_copy_n_async(scores_device.data(), numberOfAlignments, hostBuffers->scores.data(), stream);
        gw::cudautils::device_copy_n_async(pacbio_cigars_device.data(), cigarBufSize, hostBuffers->cigarOpsBuffer.data(), stream);

        GW_CU_CHECK_ERR(cudaStreamSynchronize(stream));
        int32_t i = completedAlignments;
        for(; i < numberOfAlignments; ++i)
        {
            const int32_t idx = hostBuffers->index[i] >> 1;
            const bool is_optimal = hostBuffers->index[i] & 1;
            const int64_t startPos = hostBuffers->starts[i] - hostBuffers->starts[completedAlignments];
            const int64_t len = hostBuffers->starts[i+1] - hostBuffers->starts[i];
            if(startPos + len > cigarBufSize)
            {
                if(i == completedAlignments)
                {
                    throw std::runtime_error("Could not allocate enough device or pinned host memory!");
                }
                break;
            }
            auto& res = results[idx];
            res.cigar.clear();
            res.cigar.insert(end(res.cigar), begin(hostBuffers->cigarOpsBuffer) + startPos, begin(hostBuffers->cigarOpsBuffer) + startPos + len);
            res.lastQueryPos  = querySpans[idx];
            res.lastTargetPos = targetSpans[idx];
            res.maxQueryPos   = querySpans[idx];
            res.maxTargetPos  = targetSpans[idx];
            const int32_t numEq = hostBuffers->diffs[idx].x;
            const int32_t numX  = hostBuffers->diffs[idx].y;
            const int32_t numI  = hostBuffers->diffs[idx].z;
            const int32_t numD  = hostBuffers->diffs[idx].w;
            const int32_t querySpan  = numEq + numX + numI;
            const int32_t targetSpan = numEq + numX + numD;
            assert(!is_optimal || querySpan == 0 || querySpan == querySpans[idx]);
            assert(!is_optimal || targetSpan == 0 || targetSpan == targetSpans[idx]);
            res.valid = (querySpan == querySpans[idx] && targetSpan == targetSpans[idx] && is_optimal);
            res.score    = hostBuffers->scores[idx];
            res.maxScore = hostBuffers->scores[idx];
            res.zdropped = false;
            res.diffs.numEq = numEq;
            res.diffs.numX  = numX;
            res.diffs.numI  = numI;
            res.diffs.numD  = numD;
        }
        completedAlignments += i;
    }
}

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
                                 uint32_t deviceId, int64_t maxGPUMemoryCap)
    : alnParams_(alnParams)
    , gpuStream_(claraparabricks::genomeworks::make_cuda_stream())
    , gpuHostBuffers_(std::make_unique<AlignerBatchGPUHostBuffers>())
{
    aligner_ = claraparabricks::genomeworks::cudaaligner::create_aligner(
        claraparabricks::genomeworks::cudaaligner::AlignmentType::global_alignment, maxBandwidth,
        gpuStream_.get(), deviceId, maxGPUMemoryCap);
}

AlignerBatchGPU::~AlignerBatchGPU() = default;

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

StatusAddSequencePair AlignerBatchGPU::AddSequencePairForExtensionAlignment(const char* /*query*/,
                                                                            int32_t /*queryLen*/,
                                                                            const char* /*target*/,
                                                                            int32_t /*targetLen*/)
{
    std::ostringstream oss;
    oss << "The GenomeWorks Cudaaligner does not support extension alignment at this "
           "point.";
    throw std::runtime_error(oss.str());
    return {};
}

StatusAddSequencePair AlignerBatchGPU::AddSequencePairForGlobalAlignment(const char* query,
                                                                         int32_t queryLen,
                                                                         const char* target,
                                                                         int32_t targetLen)
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
    RetrieveResultsAsPacBioCigar(gpuHostBuffers_.get(), querySpans_, targetSpans_, aligner_.get(),
                                 alnResults_, alnParams_.matchScore, alnParams_.mismatchPenalty,
                                 alnParams_.gapOpen1, alnParams_.gapExtend1);
}

}  // namespace Pancake
}  // namespace PacBio
