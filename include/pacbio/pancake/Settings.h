// Authors: Ivan Sovic

#ifndef PANCAKE_SETTINGS_H
#define PANCAKE_SETTINGS_H

#include <string>

#include <pbcopper/cli2/CLI.h>

namespace PacBio {
namespace Pancake {

struct Settings
{
    int X;
    bool Y;
    std::string Z;
    size_t NumThreads;
    std::string InputFile;
    std::string OutputFile;

    Settings();
    Settings(const PacBio::CLI_v2::Results& options);
    static PacBio::CLI_v2::Interface CreateCLI();
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_SETTINGS_H