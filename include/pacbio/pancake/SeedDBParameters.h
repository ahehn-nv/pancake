// Author: Ivan Sovic

#ifndef PANCAKE_SEEDDB_PARAMETERS_H
#define PANCAKE_SEEDDB_PARAMETERS_H

#include <cstdint>

namespace PacBio {
namespace Pancake {
namespace SeedDB {

class SeedDBParameters
{
public:
    int32_t KmerSize = 30;
    int32_t MinimizerWindow = 80;
    int32_t Spacing = 0;
    bool UseHPC = false;
    int32_t MaxHPCLen = 10;
    bool UseRC = true;

    SeedDBParameters() = default;
    ~SeedDBParameters() = default;
};

}  // namespace SeedDB
}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_SEEDDB_PARAMETERS_H