// Authors: Ivan Sovic

#ifndef PANCAKE_SEQDB_INFO_SETTINGS_H
#define PANCAKE_SEQDB_INFO_SETTINGS_H

#include <pacbio/util/CommonTypes.h>
#include <cstdint>
#include <string>

#include <pbcopper/cli2/CLI.h>

namespace PacBio {
namespace Pancake {

enum class GenomicUnit
{
    bp,
    kbp,
    Mbp,
    Gbp
};

// clang-format off
struct SeqDBInfoSettings
{
    struct Defaults
    {
        static const bool HumanReadableOutput = false;
        static const GenomicUnit Unit = GenomicUnit::bp;
    };

    std::string InputSeqDB;
    bool HumanReadableOutput = Defaults::HumanReadableOutput;
    GenomicUnit Unit = Defaults::Unit;

    SeqDBInfoSettings();
    SeqDBInfoSettings(const PacBio::CLI_v2::Results& options);
    static PacBio::CLI_v2::Interface CreateCLI();
};
// clang-format on

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_SEQDB_INFO_SETTINGS_H
