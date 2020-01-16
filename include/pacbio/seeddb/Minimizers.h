// Author: Ivan Sovic

#ifndef PANCAKE_SEEDDB_MINIMIZERS_H
#define PANCAKE_SEEDDB_MINIMIZERS_H

#include <cstdint>
#include <vector>

namespace PacBio {
namespace Pancake {
namespace SeedDB {

int GenerateMinimizers(std::vector<__int128>& minimizers, const uint8_t* seq, int32_t seqLen,
                       int32_t seqOffset, int32_t seqId, int32_t k, int32_t w, bool use_rc,
                       bool useHPC, int32_t maxHPCLen);

}  // namespace SeedDB
}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_SEEDDB_MINIMIZERS_H