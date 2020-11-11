// Author: Ivan Sovic

#ifndef PANCAKE_SEEDDB_PARAMETERS_H
#define PANCAKE_SEEDDB_PARAMETERS_H

#include <cstdint>

namespace PacBio {
namespace Pancake {
namespace SeedDB {

// clang-format off
class SeedDBParameters
{
public:
    int32_t KmerSize = 30;
    int32_t MinimizerWindow = 80;
    int32_t Spacing = 0;
    bool UseHPC = false;                // This causes the input sequences from the DB to be HP-compressed.
    bool UseHPCForSeedsOnly = false;    // This takes the uncompressed sequences, and just skips HP bases when computing seeds.
    int32_t MaxHPCLen = 10;
    bool UseRC = true;

    SeedDBParameters() = default;
    ~SeedDBParameters() = default;
};
// clang-format on

}  // namespace SeedDB
}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_SEEDDB_PARAMETERS_H