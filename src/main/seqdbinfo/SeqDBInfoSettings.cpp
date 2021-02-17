// Author: Ivan Sovic

#include "SeqDBInfoSettings.h"

#include <limits>

#include <pacbio/Version.h>

namespace PacBio {
namespace Pancake {
namespace OptionNames {

// clang-format off
const CLI_v2::PositionalArgument InputSeqDB {
R"({
    "name" : "seqdb",
    "description" : "Input SeqDB."
})"};
// clang-format on

}  // namespace OptionNames

SeqDBInfoSettings::SeqDBInfoSettings() = default;

SeqDBInfoSettings::SeqDBInfoSettings(const PacBio::CLI_v2::Results& options)
    : InputSeqDB{options[OptionNames::InputSeqDB]}
{
}

PacBio::CLI_v2::Interface SeqDBInfoSettings::CreateCLI()
{
    PacBio::CLI_v2::Interface i{"pancake seqdb-info",
                                "Computes the sequence statistics in the SeqDB.",
                                PacBio::Pancake::PancakeFormattedVersion()};

    i.DisableNumThreadsOption();

    // clang-format off
    // i.AddOptionGroup("Algorithm Options", {
    // });
    i.AddPositionalArguments({
        OptionNames::InputSeqDB,
    });

    // clang-format on
    return i;
}
}  // namespace Pancake
}  // namespace PacBio
