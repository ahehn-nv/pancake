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

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_COMMON_H