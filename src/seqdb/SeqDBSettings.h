// Authors: Ivan Sovic

#ifndef PANCAKE_SEQDB_SETTINGS_H
#define PANCAKE_SEQDB_SETTINGS_H

#include <string>
#include <cstdint>

#include <pbcopper/cli2/CLI.h>

namespace PacBio {
namespace Pancake {

struct SeqDBSettings
{
    int64_t BlockSize;
    size_t NumThreads;
    std::string InputFile;
    std::string OutputFile;

    SeqDBSettings();
    SeqDBSettings(const PacBio::CLI_v2::Results& options);
    static PacBio::CLI_v2::Interface CreateCLI();
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_SEQDB_SETTINGS_H