// Author: Armin Töpfer

#ifndef PANCAKE_VERSION_H
#define PANCAKE_VERSION_H

#include <string>

namespace PacBio {
namespace Pancake {

std::string PancakeGitSha1();
std::string PancakeVersion();

inline std::string PancakeFormattedVersion()
{
    return PancakeVersion() + " (commit " + PancakeGitSha1() + ")";
}

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_VERSION_H