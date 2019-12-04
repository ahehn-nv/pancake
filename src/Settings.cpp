// Author: Derek Barnett and Ivan Sovic

#include <pacbio/pancake/Settings.h>

#include <pacbio/pancake/Version.h>

namespace PacBio {
namespace Pancake {
namespace OptionNames {

// clang-format off

const CLI_v2::Option X{
R"({
    "names" : ["x"],
    "description" : "Dummy int opt",
    "type" : "int",
    "default" : 20
})"};

const CLI_v2::Option Y{
R"({
    "names" : ["y", "y-long"],
    "description" : "Dummy bool opt.",
    "type" : "bool"
})"};

const CLI_v2::Option Z{
R"({
    "names" : ["z"],
    "description" : "Dummy string opt.",
    "type" : "string",
    "default" : "hello"
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

Settings::Settings()
{
    // init with defaults, then modify only as needed
    // makes testing components much easier
}

Settings::Settings(const PacBio::CLI_v2::Results& options)
    : X{options[OptionNames::X]}
    , Y{options[OptionNames::Y]}
    , Z(options[OptionNames::Z])
    , NumThreads{options.NumThreads()}
    , InputFile{options[OptionNames::Input]}
    , OutputFile{options[OptionNames::Output]}
{
}

PacBio::CLI_v2::Interface Settings::CreateCLI()
{
    PacBio::CLI_v2::Interface i{"pancake", "HiFi overlapping.",
                                PacBio::Pancake::PancakeFormattedVersion()};

    // clang-format off
    i.AddOptionGroup("Algorithm Options", {
        OptionNames::X,
        OptionNames::Y,
        OptionNames::Z
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