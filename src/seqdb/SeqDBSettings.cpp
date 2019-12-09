// Author: Ivan Sovic

#include "seqdb/SeqDBSettings.h"
#include <pacbio/Version.h>

namespace PacBio {
namespace Pancake {
namespace OptionNames {

// clang-format off

const CLI_v2::PositionalArgument OutputPrefix {
R"({
    "name" : "prefix",
    "description" : "The prefix of the DB files."
})"};

const CLI_v2::PositionalArgument Input {
R"({
    "name" : "<input.fasta> [...]",
    "description" : "One or more input sequence files, in FASTA or FASTQ formats."
})"};

const CLI_v2::Option IsFofn {
R"({
    "names" : ["fofn"],
    "description" : "Input is a FOFN file."
})", SeqDBSettings::Defaults::IsFofn};

const CLI_v2::Option CompressionLevel {
R"({
    "names" : ["c", "compression"],
    "description" : "Compression level for output sequences.",
    "type" : "int"
})", SeqDBSettings::Defaults::CompressionLevel};

const CLI_v2::Option BlockSize{
R"({
    "names" : ["block-size"],
    "description" : "Block size in MB. Value 0 means all sequences will be in one block.",
    "type" : "int"
})", SeqDBSettings::Defaults::BlockSize};

// clang-format on

}  // namespace OptionNames

SeqDBSettings::SeqDBSettings() = default;

SeqDBSettings::SeqDBSettings(const PacBio::CLI_v2::Results& options)
    : OutputPrefix{options[OptionNames::OutputPrefix]}
    , InputFiles{options[OptionNames::Input]}
    , IsFofn{options[OptionNames::IsFofn]}
    , NumThreads{options.NumThreads()}
    , CompressionLevel{options[OptionNames::CompressionLevel]}
    , BlockSize{options[OptionNames::BlockSize]}
{
    // Allow multiple positional input arguments.
    const auto& files = options.PositionalArguments();
    if (files.size() < 2)
        throw std::runtime_error{"Not enough input files specified, at least one required."};
    OutputPrefix = files[0];
    InputFiles.clear();
    for (size_t i = 1; i < files.size(); ++i)
        InputFiles.push_back(files[i]);
}

PacBio::CLI_v2::Interface SeqDBSettings::CreateCLI()
{
    PacBio::CLI_v2::Interface i{"pancake", "HiFi overlapping.",
                                PacBio::Pancake::PancakeFormattedVersion()};

    // clang-format off
    i.AddOptionGroup("Algorithm Options", {
        OptionNames::IsFofn,
        OptionNames::CompressionLevel,
        OptionNames::BlockSize,
    });
    i.AddPositionalArguments({
        OptionNames::OutputPrefix,
        OptionNames::Input,
    });

    // clang-format on
    return i;
}
}  // namespace Pancake
}  // namespace PacBio