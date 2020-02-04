// Author: Ivan Sovic

#include <pacbio/Version.h>
#include <pacbio/overlaphifi/OverlapHifiSettings.h>

namespace PacBio {
namespace Pancake {
namespace OptionNames {

// clang-format off

const CLI_v2::PositionalArgument TargetDBPrefix {
R"({
    "name" : "target_prefix",
    "description" : "Prefix of the target SeqDB and SeedDB files. It should match."
})"};

const CLI_v2::PositionalArgument QueryDBPrefix {
R"({
    "name" : "query_prefix",
    "description" : "Prefix of the query SeqDB and SeedDB files. It should match."
})"};

const CLI_v2::PositionalArgument TargetBlockId {
R"({
    "name" : "target_block",
    "type" : "int",
    "description" : "Block ID from the target DB. Queries will be mapped only onto this block."
})"};

const CLI_v2::PositionalArgument QueryBlockStartId {
R"({
    "name" : "query_block_start",
    "type" : "int",
    "description" : "Start block ID for a range of blocks to map. Zero based."
})"};

const CLI_v2::PositionalArgument QueryBlockEndId {
R"({
    "name" : "query_block_end",
    "type" : "int",
    "description" : "Start block ID for a range of blocks to map. Zero based, non-inclusive."
})"};



const CLI_v2::Option FreqPercentile{
R"({
    "names" : ["freq-percentile"],
    "description" : "Filter frequent kmers.",
    "type" : "double"
})", OverlapHifiSettings::Defaults::FreqPercentile};

const CLI_v2::Option MinQueryLen{
R"({
    "names" : ["min-qlen"],
    "description" : "Ignore queries shorter than this.",
    "type" : "int"
})", OverlapHifiSettings::Defaults::MinQueryLen};

const CLI_v2::Option MaxSeedDistance{
R"({
    "names" : ["max-seed-dist"],
    "description" : "Maximum distance between two seeds to join into an anchor.",
    "type" : "int"
})", OverlapHifiSettings::Defaults::MaxSeedDistance};

const CLI_v2::Option MinNumSeeds{
R"({
    "names" : ["min-num-seeds"],
    "description" : "Minimum number of seeds in an anchor.",
    "type" : "int"
})", OverlapHifiSettings::Defaults::MinNumSeeds};

const CLI_v2::Option MinCoveredBases{
R"({
    "names" : ["min-cov-bases"],
    "description" : "Minimum number of bases covered by kmers in an anchor.",
    "type" : "int"
})", OverlapHifiSettings::Defaults::MinCoveredBases};

const CLI_v2::Option MinChainSpan{
R"({
    "names" : ["min-anchor-span"],
    "description" : "Minimum chain span to retain it.",
    "type" : "int"
})", OverlapHifiSettings::Defaults::MinChainSpan};

const CLI_v2::Option ChainBandwidth{
R"({
    "names" : ["chain-bw"],
    "description" : "Diagonal bandwidth to merge seeds into chains.",
    "type" : "int"
})", OverlapHifiSettings::Defaults::ChainBandwidth};

const CLI_v2::Option AlignmentBandwidth{
R"({
    "names" : ["aln-bw"],
    "description" : "Bandwidth for alignment, fraction of the query span.",
    "type" : "double"
})", OverlapHifiSettings::Defaults::AlignmentBandwidth};

const CLI_v2::Option AlignmentMaxD{
R"({
    "names" : ["aln-diff-rate"],
    "description" : "Expected maximum diff rate between sequences.",
    "type" : "double"
})", OverlapHifiSettings::Defaults::AlignmentMaxD};

const CLI_v2::Option MinIdentity{
R"({
    "names" : ["min-idt"],
    "description" : "Minimum percent alignment identity allowed to report the alignment.",
    "type" : "double"
})", OverlapHifiSettings::Defaults::MinIdentity};

const CLI_v2::Option MinMappedLength{
R"({
    "names" : ["min-map-len"],
    "description" : "Output only alignments above this length.",
    "type" : "int"
})", OverlapHifiSettings::Defaults::MinMappedLength};

const CLI_v2::Option AddSymmetricOverlaps{
R"({
    "names" : ["add-sym-arcs"],
    "description" : "For every overlap, it's reverse complement is also written.",
    "type" : "bool"
})", OverlapHifiSettings::Defaults::AddSymmetricOverlaps};

const CLI_v2::Option OneHitPerTarget{
R"({
    "names" : ["one-hit-per-target"],
    "description" : "Allow only one alignment per query/target pair.",
    "type" : "bool"
})", OverlapHifiSettings::Defaults::OneHitPerTarget};

// clang-format on

}  // namespace OptionNames

OverlapHifiSettings::OverlapHifiSettings() = default;

OverlapHifiSettings::OverlapHifiSettings(const PacBio::CLI_v2::Results& options)
    : TargetDBPrefix{options[OptionNames::TargetDBPrefix]}
    , QueryDBPrefix{options[OptionNames::QueryDBPrefix]}
    , NumThreads{options.NumThreads()}
    , TargetBlockId{std::stoi(options[OptionNames::TargetBlockId])}
    , QueryBlockStartId{std::stoi(options[OptionNames::QueryBlockStartId])}
    , QueryBlockEndId{std::stoi(options[OptionNames::QueryBlockEndId])}

    , FreqPercentile{options[OptionNames::FreqPercentile]}
    , MinQueryLen{options[OptionNames::MinQueryLen]}
    , MaxSeedDistance{options[OptionNames::MaxSeedDistance]}
    , MinNumSeeds{options[OptionNames::MinNumSeeds]}
    , MinCoveredBases{options[OptionNames::MinCoveredBases]}
    , MinChainSpan{options[OptionNames::MinChainSpan]}
    , ChainBandwidth{options[OptionNames::ChainBandwidth]}
    , AlignmentBandwidth{options[OptionNames::AlignmentBandwidth]}
    , AlignmentMaxD{options[OptionNames::AlignmentMaxD]}
    , MinIdentity{options[OptionNames::MinIdentity]}
    , MinMappedLength{options[OptionNames::MinMappedLength]}
    , AddSymmetricOverlaps{options[OptionNames::AddSymmetricOverlaps]}
    , OneHitPerTarget{options[OptionNames::OneHitPerTarget]}
{
}

PacBio::CLI_v2::Interface OverlapHifiSettings::CreateCLI()
{
    PacBio::CLI_v2::Interface i{"pancake", "HiFi overlapping.",
                                PacBio::Pancake::PancakeFormattedVersion()};

    // clang-format off
    i.AddOptionGroup("Algorithm Options", {
        OptionNames::FreqPercentile,
        OptionNames::MinQueryLen,
        OptionNames::MaxSeedDistance,
        OptionNames::MinNumSeeds,
        OptionNames::MinCoveredBases,
        OptionNames::MinChainSpan,
        OptionNames::ChainBandwidth,
        OptionNames::AlignmentBandwidth,
        OptionNames::AlignmentMaxD,
        OptionNames::MinIdentity,
        OptionNames::MinMappedLength,
        OptionNames::AddSymmetricOverlaps,
        OptionNames::OneHitPerTarget
    });
    i.AddPositionalArguments({
        OptionNames::TargetDBPrefix,
        OptionNames::QueryDBPrefix,
        OptionNames::TargetBlockId,
        OptionNames::QueryBlockStartId,
        OptionNames::QueryBlockEndId
    });

    // clang-format on
    return i;
}
}  // namespace Pancake
}  // namespace PacBio