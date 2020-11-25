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

class MapperBatch
{
public:
    MapperBatch(const MapperCLRSettings& settings, int32_t numThreads);
    ~MapperBatch();

    std::vector<MapperBaseResult> MapAndAlign(const std::vector<FastaSequenceCached>& targetSeqs,
                                              const std::vector<FastaSequenceCached>& querySeqs);

private:
    MapperCLRSettings settings_;
    int32_t numThreads_;
    std::unique_ptr<MapperBase> mapper_;

    static std::vector<MapperBaseResult> MapAndAlignImpl_(
        const std::vector<FastaSequenceCached>& targetSeqs,
        const std::vector<FastaSequenceCached>& querySeqs,
        const std::unique_ptr<MapperBase>& mapper, MapperCLRSettings settings, int32_t numThreads);
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_MAPPER_BATCH_H
