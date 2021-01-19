// Authors: Ivan Sovic

#include <pacbio/pancake/MapperBatch.h>
#include <pacbio/util/TicToc.h>
#include <pbcopper/logging/Logging.h>
#include <algorithm>
#include <array>
#include <iostream>
#include <sstream>
#include <tuple>

namespace PacBio {
namespace Pancake {

MapperBatch::MapperBatch(const MapperCLRSettings& settings, int32_t numThreads)
    : settings_{settings}, numThreads_(numThreads), mapper_(nullptr)
{
    // Deactivate alignment for the mapper.
    MapperCLRSettings settingsCopy = settings;
    settingsCopy.align = false;
    mapper_ = std::make_unique<MapperCLR>(settingsCopy);
}

MapperBatch::~MapperBatch() = default;

std::vector<std::vector<MapperBaseResult>> MapperBatch::DummyMapAndAlign(
    const std::vector<MapperBatchChunk>& batchData)
{
    return DummyMapAndAlignImpl_(mapper_, batchData, settings_, numThreads_);
}

std::vector<std::vector<MapperBaseResult>> MapperBatch::DummyMapAndAlignImpl_(
    std::unique_ptr<MapperCLR>& mapper, const std::vector<MapperBatchChunk>& batchData,
    MapperCLRSettings settings, int32_t /*numThreads*/)
{
    std::vector<std::vector<MapperBaseResult>> results;
    results.reserve(batchData.size());

    for (size_t i = 0; i < batchData.size(); ++i) {
        const auto& bd = batchData[i];
        std::vector<MapperBaseResult> result = mapper->MapAndAlign(bd.targetSeqs, bd.querySeqs);
        results.emplace_back(std::move(result));
    }

    if (settings.align) {
        for (size_t i = 0; i < results.size(); ++i) {
            auto& result = results[i];
            const auto& bd = batchData[i];
            for (size_t qId = 0; qId < result.size(); ++qId) {
                result[qId] = mapper->Align(bd.targetSeqs, bd.querySeqs[qId], result[qId]);
            }
        }
    }

    return results;
}

}  // namespace Pancake
}  // namespace PacBio
