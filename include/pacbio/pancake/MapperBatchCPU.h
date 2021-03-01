// Author: Ivan Sovic

#ifndef PANCAKE_MAPPER_BATCH_H
#define PANCAKE_MAPPER_BATCH_H

#include <pacbio/pancake/AlignerBatchCPU.h>
#include <pacbio/pancake/FastaSequenceCached.h>
#include <pacbio/pancake/MapperBase.h>
#include <pacbio/pancake/MapperBatchBase.h>
#include <pacbio/pancake/MapperBatchUtility.h>
#include <pacbio/pancake/MapperCLR.h>
#include <pbcopper/parallel/FireAndForget.h>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace PacBio {
namespace Pancake {

class MapperBatchCPU : public MapperBatchBase
{
public:
    MapperBatchCPU(const MapperCLRSettings& settings, int32_t numThreads);
    MapperBatchCPU(const MapperCLRSettings& settings, Parallel::FireAndForget* faf);

    ~MapperBatchCPU();

    std::vector<std::vector<MapperBaseResult>> DummyMapAndAlign(
        const std::vector<MapperBatchChunk>& batchData);

    std::vector<std::vector<MapperBaseResult>> MapAndAlign(
        const std::vector<MapperBatchChunk>& batchData) override;

private:
    MapperCLRSettings settings_;
    Parallel::FireAndForget* faf_;
    std::unique_ptr<Parallel::FireAndForget> fafFallback_;

    static std::vector<std::vector<MapperBaseResult>> DummyMapAndAlignImpl_(
        const std::vector<MapperBatchChunk>& batchChunks, MapperCLRSettings settings,
        Parallel::FireAndForget* faf);

    static std::vector<std::vector<MapperBaseResult>> MapAndAlignImpl_(
        const std::vector<MapperBatchChunk>& batchChunks, MapperCLRSettings settings,
        Parallel::FireAndForget* faf);

    static void WorkerMapper_(const std::vector<MapperBatchChunk>& batchChunks, int32_t startId,
                              int32_t endId, std::vector<std::vector<MapperBaseResult>>& results);
};

void UpdateSecondaryAndFilter(std::vector<std::vector<MapperBaseResult>>& mappingResults,
                              double secondaryAllowedOverlapFractionQuery,
                              double secondaryAllowedOverlapFractionTarget,
                              double secondaryMinScoreFraction, int32_t bestNSecondary,
                              Parallel::FireAndForget* faf);

void UpdateSecondaryAndFilter(std::vector<std::vector<MapperBaseResult>>& mappingResults,
                              Parallel::FireAndForget* faf,
                              const std::vector<MapperBatchChunk>& batchChunks);

int32_t AlignPartsOnCpu(const AlignerType& alignerTypeGlobal,
                        const AlignmentParameters& alnParamsGlobal,
                        const AlignerType& alignerTypeExt, const AlignmentParameters& alnParamsExt,
                        const std::vector<PairForBatchAlignment>& parts, const int32_t numThreads,
                        std::vector<AlignmentResult>& retAlns);

int32_t AlignPartsOnCpu(const AlignerType& alignerTypeGlobal,
                        const AlignmentParameters& alnParamsGlobal,
                        const AlignerType& alignerTypeExt, const AlignmentParameters& alnParamsExt,
                        const std::vector<PairForBatchAlignment>& parts,
                        Parallel::FireAndForget* faf, std::vector<AlignmentResult>& retAlns);

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_MAPPER_BATCH_H
