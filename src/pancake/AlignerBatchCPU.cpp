// Authors: Ivan Sovic

#include <pacbio/pancake/AlignerBatchCPU.h>
#include <pacbio/util/Util.h>
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
{
}

AlignerBatchCPU::~AlignerBatchCPU() {}

void AlignerBatchCPU::Clear()
{
    queries_.clear();
    targets_.clear();
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
    queries_.emplace_back(std::string(query, queryLen));
    targets_.emplace_back(std::string(target, targetLen));
    isGlobalAlignment_.emplace_back(isGlobalAlignment);
    return StatusAddSequencePair::OK;
}

void AlignerBatchCPU::AlignAll(int32_t numThreads)
{
    if (queries_.size() != targets_.size()) {
        std::ostringstream oss;
        oss << "Number of query and target ranges does not match. queries_.size() = "
            << queries_.size() << ", targets_.size() = " << targets_.size();
        throw std::runtime_error(oss.str());
    }
    if (queries_.size() != isGlobalAlignment_.size()) {
        std::ostringstream oss;
        oss << "Number of isGlobalAlignment_ elements does not match the number of ranges. "
               "queries_.size() = "
            << queries_.size() << ", isGlobalAlignment_.size() = " << isGlobalAlignment_.size();
        throw std::runtime_error(oss.str());
    }

    const int32_t numAlns = queries_.size();

    alnResults_.clear();
    alnResults_.resize(numAlns);

    // Nothing to do.
    if (numAlns == 0) {
        return;
    }

    // Determine how many records should land in each thread, spread roughly evenly.
    const int32_t numRecords = numAlns;
    const std::vector<std::pair<int32_t, int32_t>> jobsPerThread =
        PacBio::Pancake::DistributeJobLoad<int32_t>(numThreads, numRecords);

    // Initialize the aligners for each thread.
    std::vector<AlignerBasePtr> alignersGlobal;
    std::vector<AlignerBasePtr> alignersExt;
    for (size_t i = 0; i < jobsPerThread.size(); ++i) {
        AlignerBasePtr alignerGlobal = AlignerFactory(alnTypeGlobal_, alnParamsGlobal_);
        alignersGlobal.emplace_back(alignerGlobal);
        AlignerBasePtr alignerExt = AlignerFactory(alnTypeExt_, alnParamsExt_);
        alignersExt.emplace_back(alignerExt);
    }

    // Run the mapping in parallel.
    PacBio::Parallel::FireAndForget faf(numThreads);
    for (size_t i = 0; i < jobsPerThread.size(); ++i) {
        const int32_t jobStart = jobsPerThread[i].first;
        const int32_t jobEnd = jobsPerThread[i].second;
        faf.ProduceWith(Worker_, std::cref(queries_), std::cref(targets_),
                        std::cref(isGlobalAlignment_), jobStart, jobEnd,
                        std::ref(alignersGlobal[i]), std::ref(alignersExt[i]),
                        std::ref(alnResults_));
    }
    faf.Finalize();
}

void AlignerBatchCPU::Worker_(const std::vector<std::string>& queries,
                              const std::vector<std::string>& targets,
                              const std::vector<bool>& isGlobalAlignment, int32_t alnStartId,
                              int32_t alnEndId, AlignerBasePtr& alignerGlobal,
                              AlignerBasePtr& alignerExt, std::vector<AlignmentResult>& alnResults)
{

    for (int32_t alnId = alnStartId; alnId < alnEndId; ++alnId) {
        const auto& query = queries[alnId];
        const auto& target = targets[alnId];
        const bool isGlobal = isGlobalAlignment[alnId];

        AlignmentResult alnRes;
        if (isGlobal) {
            alnRes =
                alignerGlobal->Global(query.c_str(), query.size(), target.c_str(), target.size());
        } else {
            alnRes = alignerExt->Extend(query.c_str(), query.size(), target.c_str(), target.size());
        }
        alnResults[alnId] = std::move(alnRes);
    }
}

}  // namespace Pancake
}  // namespace PacBio
