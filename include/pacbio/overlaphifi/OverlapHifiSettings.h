// Authors: Ivan Sovic

#ifndef PANCAKE_OVERLAP_HIFI_SETTINGS_H
#define PANCAKE_OVERLAP_HIFI_SETTINGS_H

#include <cstdint>
#include <string>

#include <pbcopper/cli2/CLI.h>

namespace PacBio {
namespace Pancake {

struct OverlapHifiSettings
{
    struct Defaults
    {
        static const size_t NumThreads = 1;
        static const int32_t TargetBlockId = 0;
        static const int32_t QueryBlockStartId = 0;
        static const int32_t QueryBlockEndId = 0;

        static constexpr double FreqPercentile = 0.0002;
        static const int64_t MinQueryLen = 50;
        static const int64_t MinTargetLen = 50;
        static const int64_t MaxSeedDistance = 5000;
        static const int64_t MinNumSeeds = 3;
        static const int64_t MinCoveredBases = 30;
        static const int64_t MinChainSpan = 1000;
        static const int64_t ChainBandwidth = 100;
        static constexpr double AlignmentBandwidth = 0.01;
        static constexpr double AlignmentMaxD = 0.03;
        static constexpr double MinIdentity = 98.0;
        static const int64_t MinMappedLength = 1000;
        static const bool SkipSymmetricOverlaps = false;
        static const bool OneHitPerTarget = false;
        static const bool WriteReverseOverlaps = false;
        static const bool WriteIds = false;
        static const int32_t AllowedDovetailDist = 0;
    };

    std::string TargetDBPrefix;
    std::string QueryDBPrefix;
    size_t NumThreads = Defaults::NumThreads;

    int32_t TargetBlockId = Defaults::TargetBlockId;
    int32_t QueryBlockStartId = Defaults::QueryBlockStartId;
    int32_t QueryBlockEndId = Defaults::QueryBlockEndId;

    double FreqPercentile = Defaults::FreqPercentile;
    int64_t MinQueryLen = Defaults::MinQueryLen;
    int64_t MinTargetLen = Defaults::MinTargetLen;
    int64_t MaxSeedDistance = Defaults::MaxSeedDistance;
    int64_t MinNumSeeds = Defaults::MinNumSeeds;
    int64_t MinCoveredBases = Defaults::MinCoveredBases;
    int64_t MinChainSpan = Defaults::MinChainSpan;
    int64_t ChainBandwidth = Defaults::ChainBandwidth;
    double AlignmentBandwidth = Defaults::AlignmentBandwidth;
    double AlignmentMaxD = Defaults::AlignmentMaxD;
    double MinIdentity = Defaults::MinIdentity;
    int64_t MinMappedLength = Defaults::MinMappedLength;
    bool SkipSymmetricOverlaps = Defaults::SkipSymmetricOverlaps;
    bool OneHitPerTarget = Defaults::OneHitPerTarget;
    bool WriteReverseOverlaps = Defaults::WriteReverseOverlaps;
    bool WriteIds = Defaults::WriteIds;
    int32_t AllowedDovetailDist = Defaults::AllowedDovetailDist;

    OverlapHifiSettings();
    OverlapHifiSettings(const PacBio::CLI_v2::Results& options);
    static PacBio::CLI_v2::Interface CreateCLI();
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_SEQDB_SETTINGS_H
