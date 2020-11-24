// Author: Ivan Sovic

#ifndef PANCAKE_MAPPER_BASE_H
#define PANCAKE_MAPPER_BASE_H

#include <pacbio/pancake/DPChain.h>
#include <pacbio/pancake/FastaSequenceCached.h>
#include <pacbio/pancake/FastaSequenceId.h>
#include <pacbio/pancake/Overlap.h>
#include <pacbio/pancake/SeedIndex.h>
#include <cstdint>
#include <memory>
#include <unordered_map>
#include <vector>

namespace PacBio {
namespace Pancake {

struct ChainedRegion
{
    ChainedHits chain;
    OverlapPtr mapping;
    // Priority 0 means primary alignment, 1 secondary, < 0 not set, and > 1 filtered.
    int32_t priority = 0;
    bool isSupplementary = false;
};

class MapperBaseResult
{
public:
    std::vector<std::unique_ptr<ChainedRegion>> mappings;
};

class MapperBase
{
public:
    virtual ~MapperBase() = default;

    virtual std::vector<MapperBaseResult> MapAndAlign(
        const std::vector<std::string>& targetSeqs, const std::vector<std::string>& querySeqs) = 0;

    virtual std::vector<MapperBaseResult> MapAndAlign(
        const std::vector<FastaSequenceId>& targetSeqs,
        const std::vector<FastaSequenceId>& querySeqs) = 0;

    virtual std::vector<MapperBaseResult> MapAndAlign(
        const std::vector<FastaSequenceCached>& targetSeqs,
        const std::vector<FastaSequenceCached>& querySeqs) = 0;

    virtual MapperBaseResult MapAndAlign(const std::vector<FastaSequenceCached>& targetSeqs,
                                         const PacBio::Pancake::SeedIndex& index,
                                         const FastaSequenceCached& querySeq,
                                         const std::vector<PacBio::Pancake::Int128t>& querySeeds,
                                         const int32_t queryId, int64_t freqCutoff) = 0;
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_MAPPER_CLR_H
