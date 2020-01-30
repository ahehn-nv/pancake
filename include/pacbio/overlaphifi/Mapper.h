// Author: Ivan Sovic

#ifndef PANCAKE_OVERLAPHIFI_OVERLAPPER_H
#define PANCAKE_OVERLAPHIFI_OVERLAPPER_H

#include <pacbio/overlaphifi/OverlapHifiSettings.h>
#include <pacbio/overlaphifi/SeedIndex.h>
#include <pacbio/seeddb/Seed.h>
#include <pacbio/seeddb/SeedDBIndexCache.h>
#include <pacbio/seeddb/SequenceSeeds.h>
#include <pacbio/seqdb/FastaSequenceId.h>
#include <pacbio/seqdb/SeqDBReaderCached.h>
#include <cstdint>
#include <memory>
#include <unordered_map>
#include <vector>

namespace PacBio {
namespace Pancake {

class Mapper
{
public:
    Mapper(const OverlapHifiSettings& settings) : settings_(settings) {}
    ~Mapper() = default;

    void Map(const PacBio::Pancake::SeqDBReaderCached& targetSeqs,
             const PacBio::Pancake::SeedIndex& index,
             const PacBio::Pancake::FastaSequenceId& querySeq,
             const PacBio::Pancake::SequenceSeeds& querySeeds, int64_t freqCutoff) const;

private:
    OverlapHifiSettings settings_;

    static void DebugWriteSeedHits_(const std::string& outPath, const std::vector<SeedHit>& hits,
                                    int32_t seedLen, const std::string& queryName, int64_t queryLen,
                                    const std::string& targetName, int64_t targetLen);
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_SEQDB_READER_H