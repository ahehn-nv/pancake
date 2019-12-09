// Author: Ivan Sovic

#include "seqdb/SeqDBSettings.h"
#include <pacbio/pancake/Version.h>

namespace PacBio {
namespace Pancake {
namespace OptionNames {

// clang-format off

const CLI_v2::Option BlockSize{
R"({
    "names" : ["block-size"],
    "description" : "Block size in MB.",
    "type" : "int",
    "default" : 1000
})"};

const CLI_v2::PositionalArgument Input {
R"({
    "name" : "input.*",
    "description" : "Input data"
})"};

const CLI_v2::PositionalArgument Output {
R"({
    "name" : "output.fasta",
    "description" : "Output results"
})"};

// clang-format on

}  // namespace OptionNames

SeqDBSettings::SeqDBSettings()
{
    // init with defaults, then modify only as needed
    // makes testing components much easier
}

SeqDBSettings::SeqDBSettings(const PacBio::CLI_v2::Results& options)
    : BlockSize{options[OptionNames::BlockSize]}
    , NumThreads{options.NumThreads()}
    , InputFile{options[OptionNames::Input]}
    , OutputFile{options[OptionNames::Output]}
{
}

PacBio::CLI_v2::Interface SeqDBSettings::CreateCLI()
{
    PacBio::CLI_v2::Interface i{"pancake", "HiFi overlapping.",
                                PacBio::Pancake::PancakeFormattedVersion()};

    // clang-format off
    i.AddOptionGroup("Algorithm Options", {
        OptionNames::BlockSize
    });
    i.AddPositionalArguments({
        OptionNames::Input,
        OptionNames::Output
    });

    // clang-format on
    return i;
}
}  // namespace Pancake
}  // namespace PacBio