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

const CLI_v2::Option HumanReadableOutput {
R"({
    "names" : ["human"],
    "description" : "Write human-readable output instead of Json.",
    "type" : "bool"
})", SeqDBInfoSettings::Defaults::HumanReadableOutput};

const CLI_v2::Option GenomicUnit {
R"({
    "names" : ["unit"],
    "choices" : ["bp", "kbp", "Mbp", "Gbp"],
    "description" : "Unit for genomic length for the output stats.",
    "type" : "string"
})", "bp"};

// clang-format on

}  // namespace OptionNames

GenomicUnit GenomicUnitFromString(const std::string& val)
{
    if (val == "bp") {
        return GenomicUnit::bp;
    } else if (val == "kbp") {
        return GenomicUnit::kbp;
    } else if (val == "Mbp") {
        return GenomicUnit::Mbp;
    } else if (val == "Gbp") {
        return GenomicUnit::Gbp;
    }
    throw std::runtime_error("Unknown genomic unit: '" + val + "'.");
}

SeqDBInfoSettings::SeqDBInfoSettings() = default;

SeqDBInfoSettings::SeqDBInfoSettings(const PacBio::CLI_v2::Results& options)
    : InputSeqDB{options[OptionNames::InputSeqDB]}
    , HumanReadableOutput(options[OptionNames::HumanReadableOutput])
    , Unit(GenomicUnitFromString(options[OptionNames::GenomicUnit]))
{
}

PacBio::CLI_v2::Interface SeqDBInfoSettings::CreateCLI()
{
    PacBio::CLI_v2::Interface i{"pancake seqdb-info",
                                "Computes the sequence statistics in the SeqDB.",
                                PacBio::Pancake::PancakeFormattedVersion()};

    i.DisableNumThreadsOption();

    // clang-format off
    i.AddOptionGroup("Options", {
        OptionNames::HumanReadableOutput,
        OptionNames::GenomicUnit,
    });
    i.AddPositionalArguments({
        OptionNames::InputSeqDB,
    });

    // clang-format on
    return i;
}
}  // namespace Pancake
}  // namespace PacBio
