// Authors: Ivan Sovic

#ifndef PANCAKE_SEEDDB_SETTINGS_H
#define PANCAKE_SEEDDB_SETTINGS_H

#include <cstdint>
#include <string>

#include <pbcopper/cli2/CLI.h>

namespace PacBio {
namespace Pancake {
namespace SeedDB {

struct SeedDBSettings
{
    struct Defaults
    {
        static const size_t NumThreads = 1;
        static constexpr float BufferSize = 1000.0f;
        static const bool SplitBlocks = false;
        static const int32_t KmerSize = 30;
        static const int32_t MinimizerWindow = 80;
        static const bool UseHPC = false;
        static const int32_t MaxHPCLen = 10;
    };

    std::string InputFile;
    std::string OutputPrefix;
    size_t NumThreads = Defaults::NumThreads;
    float BufferSize = Defaults::BufferSize;
    bool SplitBlocks = Defaults::SplitBlocks;
    int32_t KmerSize = Defaults::KmerSize;
    int32_t MinimizerWindow = Defaults::MinimizerWindow;
    bool UseHPC = Defaults::UseHPC;
    int32_t MaxHPCLen = Defaults::MaxHPCLen;

    SeedDBSettings();
    SeedDBSettings(const PacBio::CLI_v2::Results& options);
    static PacBio::CLI_v2::Interface CreateCLI();
};

}  // namespace SeedDB
}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_SEQDB_SETTINGS_H