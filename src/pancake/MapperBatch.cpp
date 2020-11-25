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
    : settings_{settings}, numThreads_(numThreads), mapper_(std::make_unique<MapperCLR>(settings))
{
}

MapperBatch::~MapperBatch() = default;

std::vector<MapperBaseResult> MapperBatch::MapAndAlign(
    const std::vector<FastaSequenceCached>& targetSeqs,
    const std::vector<FastaSequenceCached>& querySeqs)
{
    MapAndAlignImpl_(targetSeqs, querySeqs, mapper_, settings_, numThreads_);
}

std::vector<MapperBaseResult> MapperBatch::MapAndAlignImpl_(
    const std::vector<FastaSequenceCached>& targetSeqs,
    const std::vector<FastaSequenceCached>& querySeqs, const std::unique_ptr<MapperBase>& mapper,
    MapperCLRSettings settings, int32_t numThreads)
{
    return {};
}

}  // namespace Pancake
}  // namespace PacBio
