// Author: Ivan Sovic

#ifndef PANCAKE_SEEDDB_MINIMIZERS_H
#define PANCAKE_SEEDDB_MINIMIZERS_H

#include <pacbio/pancake/FastaSequenceCached.h>
#include <pacbio/pancake/Seed.h>
#include <pacbio/pancake/SeedHit.h>
#include <pacbio/util/CommonTypes.h>
#include <array>
#include <cstdint>
#include <deque>
#include <vector>

namespace PacBio {
namespace Pancake {
namespace SeedDB {

static inline uint64_t InvertibleHash(uint64_t key, uint64_t mask)
{
    /*
    Credit: Heng Li, Minimap2.
    */
    key = (~key + (key << 21)) & mask;  // key = (key << 21) - key - 1;
    key = key ^ key >> 24;
    key = ((key + (key << 3)) + (key << 8)) & mask;  // key * 265
    key = key ^ key >> 14;
    key = ((key + (key << 2)) + (key << 4)) & mask;  // key * 21
    key = key ^ key >> 28;
    key = (key + (key << 31)) & mask;
    return key;
}

static inline uint64_t ComputeKmerMask(int32_t kmerSize)
{
    const uint64_t mask =
        (kmerSize < 32) ? ((((uint64_t)1) << (2 * kmerSize)) - 1)
                        : 0xFFFFFFFFFFFFFFFF;  // Mask the number of required bits for the kmer.
    return mask;
}

int GenerateMinimizers(std::vector<PacBio::Pancake::Int128t>& minimizers, const uint8_t* seq,
                       const int32_t seqLen, const int32_t seqOffset, const int32_t seqId,
                       const int32_t kmerSize, const int32_t winSize, const int32_t spacing,
                       const bool useReverseComplement, const bool useHPC, const int32_t maxHPCLen);

void GenerateMinimizers(std::vector<PacBio::Pancake::Int128t>& retSeeds,
                        std::vector<int32_t>& retSequenceLengths,
                        const std::vector<FastaSequenceCached>& targetSeqs, const int32_t kmerSize,
                        const int32_t winSize, const int32_t spacing,
                        const bool useReverseComplement, const bool useHPC,
                        const int32_t maxHPCLen);

void GenerateMinimizers(std::vector<PacBio::Pancake::Int128t>& retSeeds,
                        std::vector<int32_t>& retSequenceLengths,
                        const std::vector<std::string>& targetSeqs, const int32_t kmerSize,
                        const int32_t winSize, const int32_t spacing,
                        const bool useReverseComplement, const bool useHPC,
                        const int32_t maxHPCLen);

template <class TargetHashType>
bool CollectSeedHits(std::vector<SeedHit>& hits, const PacBio::Pancake::SeedDB::SeedRaw* querySeeds,
                     const int64_t querySeedsSize, const int32_t /*queryLen*/,
                     const TargetHashType& hash,
                     const PacBio::Pancake::SeedDB::SeedRaw* targetSeeds,
                     const int64_t /*targetSeedsSize*/, const std::vector<int32_t>& targetLengths,
                     const int32_t /*kmerSize*/, const int32_t /*spacing*/,
                     const int64_t freqCutoff)
{
    hits.clear();

    static auto GetSequenceLength = [](const std::vector<int32_t>& sequenceLengths,
                                       int32_t seqId) -> int32_t {
        // Sanity check for the sequence ID.
        if (seqId < 0 || seqId >= static_cast<int32_t>(sequenceLengths.size())) {
            std::ostringstream oss;
            oss << "Invalid seqId. seqId = " << seqId
                << ", sequenceLengths.size() = " << sequenceLengths.size();
            throw std::runtime_error(oss.str());
        }
        return sequenceLengths[seqId];
    };

    // The +1 is because for every seed base there are Spacing spaces, and the subtraction
    // is because after the last seed base the spaces shouldn't be counted.
    //    const int32_t seedSize = kmerSize * (spacing + 1) - spacing;

    for (int64_t seedId = 0; seedId < querySeedsSize; ++seedId) {
        const auto& querySeed = querySeeds[seedId];
        auto decodedQuery = PacBio::Pancake::SeedDB::Seed(querySeed);
        auto it = hash.find(decodedQuery.key);
        if (it != hash.end()) {
            int64_t start = std::get<0>(it->second);
            int64_t end = std::get<1>(it->second);
            // Skip very frequent seeds.
            if (freqCutoff > 0 && (end - start) > freqCutoff) {
                continue;
            }
            for (int64_t i = start; i < end; ++i) {
                auto decodedTarget = PacBio::Pancake::SeedDB::Seed(targetSeeds[i]);
                bool isRev = false;
                int32_t targetPos = decodedTarget.pos;  // Start position of the target kmer hit.
                int32_t queryPos = decodedQuery.pos;    // Start position of the query kmer hit.
                const int32_t querySpan = decodedQuery.span;
                const int32_t targetSpan = decodedTarget.span;

                if (decodedQuery.IsRev() != decodedTarget.IsRev()) {
                    isRev = true;
                    const int32_t targetLen = GetSequenceLength(targetLengths, decodedTarget.seqID);
                    targetPos = targetLen - (decodedTarget.pos + targetSpan);
                    // queryPos = queryLen - (decodedQuery.pos +
                    //                        querySpan);  // End pos in fwd is start pos in rev.
                }

                SeedHit hit{decodedTarget.seqID, isRev,     targetPos, queryPos,
                            targetSpan,          querySpan, 0};

                hits.emplace_back(hit);
            }
        }
    }

    return hits.size() > 0;
}

}  // namespace SeedDB
}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_SEEDDB_MINIMIZERS_H
