// Authors: Ivan Sovic

#include <pacbio/overlaphifi/Mapper.h>
#include <pacbio/util/TicToc.h>
#include <pbcopper/logging/Logging.h>
#include <iostream>
#include <sstream>

namespace PacBio {
namespace Pancake {

// #define PANCAKE_DEBUG

void Mapper::DebugWriteSeedHits_(const std::string& outPath, const std::vector<SeedHit>& hits,
                                 int32_t seedLen, const std::string& queryName, int64_t queryLen,
                                 const std::string& targetName, int64_t targetLen)
{
    std::ofstream ofs(outPath);
    // Simply walk away if the file cannot be open.
    // Avoid writing debug output if a specific path is not available.
    if (ofs.is_open() == false) {
        return;
    }
    ofs << queryName << "\t0\t" << queryLen << "\t" << targetName << "\t0\t" << targetLen << "\t0.0"
        << std::endl;
    for (size_t i = 0; i < hits.size(); ++i) {
        int32_t clusterId = hits[i].targetId * 2 + (hits[i].targetRev ? 1 : 0);
        ofs << hits[i].queryPos << "\t" << hits[i].targetPos << "\t" << clusterId << std::endl;
        ofs << hits[i].queryPos + seedLen << "\t" << hits[i].targetPos + seedLen << "\t"
            << clusterId << std::endl;
    }
}

void Mapper::Map(const PacBio::Pancake::SeqDBReaderCached& /*targetSeqs*/,
                 const PacBio::Pancake::SeedIndex& index,
                 const PacBio::Pancake::FastaSequenceId& querySeq,
                 const PacBio::Pancake::SequenceSeeds& querySeeds, int64_t freqCutoff) const
{
#ifdef PANCAKE_DEBUG
    PBLOG_INFO << "Mapping query ID = " << querySeq.Id() << ", header = " << querySeq.Name();
#endif

    if (static_cast<int64_t>(querySeq.Bases().size()) < settings_.MinQueryLen) {
        return;
    }

    TicToc ttCollectHits;
    std::vector<SeedHit> hits;
    index.CollectHits(querySeeds.Seeds(), hits, freqCutoff);
    ttCollectHits.Stop();

    TicToc ttSortHits;
    std::sort(hits.begin(), hits.end());
    ttSortHits.Stop();

#ifdef PANCAKE_DEBUG
    PBLOG_INFO << "Collected " << hits.size() << " hits.";
    PBLOG_INFO << "Time - collecting hits: " << ttCollectHits.GetMillisecs() << " ms / "
               << ttCollectHits.GetCpuMillisecs() << " CPU ms";
    PBLOG_INFO << "Time - sorting: " << ttSortHits.GetMillisecs() << " ms / "
               << ttSortHits.GetCpuMillisecs() << " CPU ms";
    DebugWriteSeedHits_("temp/debug/mapper-0-seed_hits.csv", hits, 30, querySeq.Name(),
                        querySeq.Bases().size(), "target", 0);
#endif
}

}  // namespace Pancake
}  // namespace PacBio
