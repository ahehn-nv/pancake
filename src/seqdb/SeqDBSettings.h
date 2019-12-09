// Authors: Ivan Sovic

#ifndef PANCAKE_SEQDB_SETTINGS_H
#define PANCAKE_SEQDB_SETTINGS_H

#include <cstdint>
#include <string>

#include <pbcopper/cli2/CLI.h>

namespace PacBio {
namespace Pancake {

struct SeqDBSettings
{
    struct Defaults
    {
        static const bool IsFofn = false;
        static const size_t NumThreads = 1;
        static const int32_t CompressionLevel = 1;
        static const int64_t BlockSize = 1000;
    };

    std::string OutputPrefix;
    std::vector<std::string> InputFiles;
    bool IsFofn = Defaults::IsFofn;
    size_t NumThreads = Defaults::NumThreads;
    int32_t CompressionLevel = Defaults::CompressionLevel;
    int64_t BlockSize = Defaults::BlockSize;

    SeqDBSettings();
    SeqDBSettings(const PacBio::CLI_v2::Results& options);
    static PacBio::CLI_v2::Interface CreateCLI();
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_SEQDB_SETTINGS_H