// Author: Ivan Sovic

#ifndef PANCAKE_MAPPER_BATCH_H
#define PANCAKE_MAPPER_BATCH_H

#include <pacbio/pancake/FastaSequenceCached.h>
#include <pacbio/pancake/MapperBase.h>
#include <pacbio/pancake/MapperCLR.h>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace PacBio {
namespace Pancake {

struct MapperBatchChunk
{
    std::vector<FastaSequenceCached> targetSeqs;
    std::vector<FastaSequenceCached> querySeqs;
};

class MapperBatch
{
public:
    MapperBatch(const MapperCLRSettings& settings, int32_t numThreads);
    ~MapperBatch();

    std::vector<std::vector<MapperBaseResult>> DummyMapAndAlign(
        const std::vector<MapperBatchChunk>& batchData);

private:
    MapperCLRSettings settings_;
    int32_t numThreads_;
    std::unique_ptr<MapperCLR> mapper_;

    static std::vector<std::vector<MapperBaseResult>> DummyMapAndAlignImpl_(
        std::unique_ptr<MapperCLR>& mapper, const std::vector<MapperBatchChunk>& batchData,
        MapperCLRSettings settings, int32_t numThreads);
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_MAPPER_BATCH_H
