// Author: Ivan Sovic

#ifndef PANCAKE_COMMON_H
#define PANCAKE_COMMON_H

#include <cstdint>
#include <lib/flat_hash_map/flat_hash_map.hpp>
#include <string>

namespace PacBio {
namespace Pancake {

using HeaderLookupType = ska::flat_hash_map<std::string, int32_t>;
using IdLookupType = ska::flat_hash_map<int32_t, int32_t>;

static const __int128 MASK_U64_LOW_2BIT = 0x0000000000000003;

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_COMMON_H