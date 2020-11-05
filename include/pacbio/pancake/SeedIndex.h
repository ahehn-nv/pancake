// Author: Ivan Sovic

#ifndef PANCAKE_OVERLAPHIFI_SEEDINDEX_H
#define PANCAKE_OVERLAPHIFI_SEEDINDEX_H

#include <pacbio/pancake/Seed.h>
#include <pacbio/pancake/SeedDBIndexCache.h>
#include <pacbio/pancake/SeedHit.h>
#include <cstdint>
#include <memory>
#include <unordered_map>
#include <vector>

// struct SeedHasher {
//     size_t operator()(const uint64_t& x) const {
//         return x;
//     }
// };

/*
 * Experimental code test performance with different hash maps.
 */
// #define SEED_INDEX_USING_DENSEHASH
// #define SEED_INDEX_USING_SPARSEHASH
// #define SEED_INDEX_USING_UNORDERED_MAP
#define SEED_INDEX_USING_FLATHASHMAP

const uint64_t SEED_INDEX_EMPTY_HASH_KEY = 0xFFFFFFFFFFFFFFFF;

#ifdef SEED_INDEX_USING_UNORDERED_MAP
#include <unordered_map>
// Key: kmer hash, Value: pair of <startId, endId> in the seeds_ vector.
typedef std::unordered_map<uint64_t, std::pair<int64_t, int64_t>, std::hash<uint64_t>> SeedHashType;
#endif

#ifdef SEED_INDEX_USING_DENSEHASH
#include <sparsehash/dense_hash_map>
using google::dense_hash_map;  // namespace where class lives by default
// Key: kmer hash, Value: pair of <startId, endId> in the seeds_ vector.
typedef dense_hash_map<uint64_t, std::pair<int64_t, int64_t>, std::hash<uint64_t>> SeedHashType;
#endif

#ifdef SEED_INDEX_USING_SPARSEHASH
#include <sparsehash/sparse_hash_map>
using google::sparse_hash_map;  // namespace where class lives by default
// Key: kmer hash, Value: pair of <startId, endId> in the seeds_ vector.
typedef sparse_hash_map<uint64_t, std::pair<int64_t, int64_t>, std::hash<uint64_t>> SeedHashType;
#endif

#ifdef SEED_INDEX_USING_FLATHASHMAP
#include <lib/flat_hash_map/flat_hash_map.hpp>
// Key: kmer hash, Value: pair of <startId, endId> in the seeds_ vector.
typedef ska::flat_hash_map<uint64_t, std::pair<int64_t, int64_t>> SeedHashType;
#endif

namespace PacBio {
namespace Pancake {

class SeedIndex
{
public:
    SeedIndex(std::shared_ptr<PacBio::Pancake::SeedDBIndexCache>& seedDBCache,
              std::vector<PacBio::Pancake::SeedDB::SeedRaw>&& seeds);
    SeedIndex(const PacBio::Pancake::SeedDB::SeedDBParameters& seedParams,
              const std::vector<int32_t>& sequenceLengths,
              std::vector<PacBio::Pancake::SeedDB::SeedRaw>&& seeds);
    ~SeedIndex();

    void ComputeFrequencyStats(double percentileCutoff, int64_t& retFreqMax, double& retFreqAvg,
                               double& retFreqMedian, int64_t& retFreqCutoff) const;
    int64_t GetSeeds(uint64_t key, std::vector<PacBio::Pancake::SeedDB::SeedRaw>& seeds) const;
    bool CollectHits(const std::vector<PacBio::Pancake::SeedDB::SeedRaw>& querySeeds,
                     int32_t queryLen, std::vector<SeedHit>& hits, int64_t freqCutoff) const;
    bool CollectHits(const PacBio::Pancake::SeedDB::SeedRaw* querySeeds, int64_t querySeedsSize,
                     int32_t queryLen, std::vector<SeedHit>& hits, int64_t freqCutoff) const;

    const std::vector<int32_t> GetSequenceLengths() const { return sequenceLengths_; }

    const PacBio::Pancake::SeedDB::SeedDBParameters& GetSeedParams() const { return seedParams_; }

    int32_t GetSequenceLength(int32_t seqId) const
    {
        // Sanity check for the sequence ID.
        if (seqId < 0 || seqId >= static_cast<int32_t>(sequenceLengths_.size())) {
            std::ostringstream oss;
            oss << "Invalid seqId. seqId = " << seqId
                << ", sequenceLengths_.size() = " << sequenceLengths_.size();
            throw std::runtime_error(oss.str());
        }
        return sequenceLengths_[seqId];
    }

private:
    std::vector<PacBio::Pancake::SeedDB::SeedRaw> seeds_;
    SeedHashType hash_;
    PacBio::Pancake::SeedDB::SeedDBParameters seedParams_;
    std::vector<int32_t> sequenceLengths_;

    void BuildHash_();
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_SEQDB_READER_H
