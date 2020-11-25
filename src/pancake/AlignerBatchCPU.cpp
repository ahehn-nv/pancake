// Authors: Ivan Sovic

#include <pacbio/pancake/AlignerBatchCPU.h>
#include <pbcopper/parallel/FireAndForget.h>
#include <pbcopper/parallel/WorkQueue.h>
#include <cstring>
#include <iostream>
#include <sstream>

namespace PacBio {
namespace Pancake {

AlignerBatchCPU::AlignerBatchCPU(const AlignerType alnTypeGlobal,
                                 const AlignmentParameters& alnParamsGlobal,
                                 const AlignerType alnTypeExt,
                                 const AlignmentParameters& alnParamsExt)
    : alnTypeGlobal_(alnTypeGlobal)
    , alnParamsGlobal_(alnParamsGlobal)
    , alnTypeExt_(alnTypeExt)
    , alnParamsExt_(alnParamsExt)
    , bases_{}
    , numBases_(0)
    , seqRangesQuery_{}
    , seqRangesTarget_{}
    , isGlobalAlignment_{}
    , alnResults_{}
{
}

AlignerBatchCPU::~AlignerBatchCPU() {}

void AlignerBatchCPU::Clear()
{
    bases_.clear();
    numBases_ = 0;
    seqRangesQuery_.clear();
    seqRangesTarget_.clear();
    isGlobalAlignment_.clear();
    alnResults_.clear();
}

StatusAddSequencePair AlignerBatchCPU::AddSequencePair(const char* query, int32_t queryLen,
                                                       const char* target, int32_t targetLen,
                                                       bool isGlobalAlignment)
{
    if (queryLen < 0 || targetLen < 0) {
        return StatusAddSequencePair::SEQUENCE_LEN_BELOW_ZERO;
    }

    int64_t currSize = bases_.size();
    while ((numBases_ + queryLen + targetLen) > currSize) {
        currSize = std::max(static_cast<int64_t>(1), currSize * 2);
    }
    if (currSize > static_cast<int64_t>(bases_.size())) {
        bases_.resize(currSize);
    }

    Range rangeQuery;
    rangeQuery.start = numBases_;
    rangeQuery.end = numBases_ + queryLen;
    seqRangesQuery_.emplace_back(std::move(rangeQuery));
    memcpy(&bases_[numBases_], query, queryLen);
    numBases_ += queryLen;

    Range rangeTarget;
    rangeTarget.start = numBases_;
    rangeTarget.end = numBases_ + targetLen;
    seqRangesTarget_.emplace_back(std::move(rangeTarget));
    memcpy(&bases_[numBases_], target, targetLen);
    numBases_ += targetLen;

    isGlobalAlignment_.emplace_back(isGlobalAlignment);

    return StatusAddSequencePair::OK;
}

void AlignerBatchCPU::AlignAll(int32_t numThreads)
{
    if (seqRangesQuery_.size() != seqRangesTarget_.size()) {
        std::ostringstream oss;
        oss << "Number of query and target ranges does not match. seqRangesQuery_.size() = "
            << seqRangesQuery_.size() << ", seqRangesTarget_.size() = " << seqRangesTarget_.size();
        throw std::runtime_error(oss.str());
    }
    if (seqRangesQuery_.size() != isGlobalAlignment_.size()) {
        std::ostringstream oss;
        oss << "Number of isGlobalAlignment_ elements does not match the number of ranges. "
               "seqRangesQuery_.size() = "
            << seqRangesQuery_.size()
            << ", isGlobalAlignment_.size() = " << isGlobalAlignment_.size();
        throw std::runtime_error(oss.str());
    }

    const int32_t numAlns = seqRangesQuery_.size();

    alnResults_.clear();
    alnResults_.resize(numAlns);

    // Determine how many records should land in each thread, spread roughly evenly.
    const int32_t numRecords = numAlns;
    const int32_t actualThreadCount = std::min(numThreads, numRecords);
    const int32_t minimumRecordsPerThreads = (numRecords / actualThreadCount);
    const int32_t remainingRecords = (numRecords % actualThreadCount);
    std::vector<int32_t> recordsPerThread(actualThreadCount, minimumRecordsPerThreads);
    for (int32_t i = 0; i < remainingRecords; ++i) {
        ++recordsPerThread[i];
    }

    // Initialize the aligners for each thread.
    std::vector<AlignerBasePtr> alignersGlobal;
    std::vector<AlignerBasePtr> alignersExt;
    for (int32_t i = 0; i < actualThreadCount; ++i) {
        AlignerBasePtr alignerGlobal = AlignerFactory(alnTypeGlobal_, alnParamsGlobal_);
        alignersGlobal.emplace_back(alignerGlobal);
        AlignerBasePtr alignerExt = AlignerFactory(alnTypeExt_, alnParamsExt_);
        alignersExt.emplace_back(alignerExt);
    }

    // Run the mapping in parallel.
    PacBio::Parallel::FireAndForget faf(numThreads);
    int32_t submittedCount = 0;
    for (int32_t i = 0; i < actualThreadCount; ++i) {
        faf.ProduceWith(Worker_, std::cref(bases_), std::cref(seqRangesQuery_),
                        std::cref(seqRangesTarget_), std::cref(isGlobalAlignment_), submittedCount,
                        submittedCount + recordsPerThread[i], std::ref(alignersGlobal[i]),
                        std::ref(alignersExt[i]), std::ref(alnResults_));
        submittedCount += recordsPerThread[i];
    }
    faf.Finalize();
}

void AlignerBatchCPU::Worker_(const std::vector<char>& bases,
                              const std::vector<Range>& seqRangesQuery,
                              const std::vector<Range>& seqRangesTarget,
                              const std::vector<bool>& isGlobalAlignment, int32_t alnStartId,
                              int32_t alnEndId, AlignerBasePtr& alignerGlobal,
                              AlignerBasePtr& alignerExt, std::vector<AlignmentResult>& alnResults)
{

    for (int32_t alnId = alnStartId; alnId < alnEndId; ++alnId) {
        const auto& qRange = seqRangesQuery[alnId];
        const auto& tRange = seqRangesTarget[alnId];
        const bool isGlobal = isGlobalAlignment[alnId];

        AlignmentResult alnRes;
        if (isGlobal) {
            alnRes = alignerGlobal->Global(bases.data() + qRange.start, qRange.Span(),
                                           bases.data() + tRange.start, tRange.Span());
        } else {
            alnRes = alignerExt->Extend(bases.data() + qRange.start, qRange.Span(),
                                        bases.data() + tRange.start, tRange.Span());
        }
        alnResults[alnId] = std::move(alnRes);
    }
}

}  // namespace Pancake
}  // namespace PacBio
