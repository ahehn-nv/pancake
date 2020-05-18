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
    "description" : "Start block ID for a range of blocks to map. Zero based, non-inclusive. Value == 0 runs until the end block."
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

const CLI_v2::Option MinTargetLen{
R"({
    "names" : ["min-tlen"],
    "description" : "Ignore targets shorter than this.",
    "type" : "int"
})", OverlapHifiSettings::Defaults::MinTargetLen};

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

const CLI_v2::Option SkipSymmetricOverlaps{
R"({
    "names" : ["skip-sym"],
    "description" : "If Aid < Bid, only compute overlap Aid->Bid and skip computing overlap for Bid->Aid.",
    "type" : "bool"
})", OverlapHifiSettings::Defaults::SkipSymmetricOverlaps};

const CLI_v2::Option OneHitPerTarget{
R"({
    "names" : ["one-hit-per-target"],
    "description" : "Allow only one alignment per query/target pair.",
    "type" : "bool"
})", OverlapHifiSettings::Defaults::OneHitPerTarget};

const CLI_v2::Option WriteReverseOverlaps{
R"({
    "names" : ["write-rev"],
    "description" : "For eveery overlap, write out its reverse complement too.",
    "type" : "bool"
})", OverlapHifiSettings::Defaults::WriteReverseOverlaps};
const CLI_v2::Option WriteIds{
R"({
    "names" : ["write-ids"],
    "description" : "Output overlaps will contain numeric IDs for the A and B reads (instead of names).",
    "type" : "bool"
})", OverlapHifiSettings::Defaults::WriteIds};
const CLI_v2::Option WriteCigar{
R"({
    "names" : ["write-cigar"],
    "description" : "Write the CIGAR string if the sensitive alignment mode is applied.",
    "type" : "bool"
})", OverlapHifiSettings::Defaults::WriteCigar};

const CLI_v2::Option AllowedDovetailDist{
R"({
    "names" : ["dt-dist"],
    "description" : "Allowed distance of an overlap from the beginning of the sequences to call the overlap a dovetail.",
    "type" : "int"
})", OverlapHifiSettings::Defaults::AllowedDovetailDist};

const CLI_v2::Option AllowedHeuristicExtendDist{
R"({
    "names" : ["ext-dist"],
    "description" : "Heuristically modify the coordinats of an overlap into a dovetail overlap if are within this distance from the edges of the reads.",
    "type" : "int"
})", OverlapHifiSettings::Defaults::AllowedHeuristicExtendDist};

const CLI_v2::Option CombineBlocks{
R"({
    "names" : ["combine"],
    "description" : "Combines this many query blocks into one larger block for processing.",
    "type" : "int"
})", OverlapHifiSettings::Defaults::CombineBlocks};

const CLI_v2::Option BestN{
R"({
    "names" : ["bestn"],
    "description" : "Output only best N alignments.",
    "type" : "int"
})", OverlapHifiSettings::Defaults::BestN};

const CLI_v2::Option UseHPC{
R"({
    "names" : ["use-hpc"],
    "description" : "Enable homopolymer compression."
})", OverlapHifiSettings::Defaults::UseHPC};

const CLI_v2::Option UseTraceback{
R"({
    "names" : ["traceback"],
    "description" : "Run alignment traceback and compute mismatches.",
    "type" : "bool"
})", OverlapHifiSettings::Defaults::UseTraceback};

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
    , MinTargetLen{options[OptionNames::MinTargetLen]}
    , MaxSeedDistance{options[OptionNames::MaxSeedDistance]}
    , MinNumSeeds{options[OptionNames::MinNumSeeds]}
    , MinCoveredBases{options[OptionNames::MinCoveredBases]}
    , MinChainSpan{options[OptionNames::MinChainSpan]}
    , ChainBandwidth{options[OptionNames::ChainBandwidth]}
    , AlignmentBandwidth{options[OptionNames::AlignmentBandwidth]}
    , AlignmentMaxD{options[OptionNames::AlignmentMaxD]}
    , MinIdentity{options[OptionNames::MinIdentity]}
    , MinMappedLength{options[OptionNames::MinMappedLength]}
    , SkipSymmetricOverlaps{options[OptionNames::SkipSymmetricOverlaps]}
    , OneHitPerTarget{options[OptionNames::OneHitPerTarget]}
    , WriteReverseOverlaps{options[OptionNames::WriteReverseOverlaps]}
    , WriteIds{options[OptionNames::WriteIds]}
    , WriteCigar{options[OptionNames::WriteCigar]}
    , AllowedDovetailDist{options[OptionNames::AllowedDovetailDist]}
    , AllowedHeuristicExtendDist{options[OptionNames::AllowedHeuristicExtendDist]}
    , CombineBlocks{options[OptionNames::CombineBlocks]}
    , BestN{options[OptionNames::BestN]}
    , UseHPC{options[OptionNames::UseHPC]}
    , UseTraceback{options[OptionNames::UseTraceback]}
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
        OptionNames::MinTargetLen,
        OptionNames::MaxSeedDistance,
        OptionNames::MinNumSeeds,
        OptionNames::MinCoveredBases,
        OptionNames::MinChainSpan,
        OptionNames::ChainBandwidth,
        OptionNames::AlignmentBandwidth,
        OptionNames::AlignmentMaxD,
        OptionNames::MinIdentity,
        OptionNames::MinMappedLength,
        OptionNames::SkipSymmetricOverlaps,
        OptionNames::OneHitPerTarget,
        OptionNames::WriteReverseOverlaps,
        OptionNames::WriteIds,
        OptionNames::WriteCigar,
        OptionNames::AllowedDovetailDist,
        OptionNames::AllowedHeuristicExtendDist,
        OptionNames::CombineBlocks,
        OptionNames::BestN,
        OptionNames::UseHPC,
        OptionNames::UseTraceback,
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
