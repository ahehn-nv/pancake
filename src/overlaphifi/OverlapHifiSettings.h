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
        static const int64_t BlockSize = 1000;
        static const size_t NumThreads = 1;
    };

    int64_t BlockSize = Defaults::BlockSize;
    size_t NumThreads = Defaults::NumThreads;
    std::string InputFile;
    std::string OutputFile;

    OverlapHifiSettings();
    OverlapHifiSettings(const PacBio::CLI_v2::Results& options);
    static PacBio::CLI_v2::Interface CreateCLI();
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_SEQDB_SETTINGS_H