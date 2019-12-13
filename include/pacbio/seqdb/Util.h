// Author: Ivan Sovic

#ifndef PANCAKE_UTIL_H
#define PANCAKE_UTIL_H

#include <cstdio>
#include <memory>
#include <sstream>

namespace PacBio {
namespace Pancake {

struct FileDeleter
{
    void operator()(std::FILE* fp) const { fclose(fp); }
};

inline std::unique_ptr<FILE, FileDeleter> OpenFile(const std::string& filename,
                                                   const std::string& mode)
{
    auto ret = std::unique_ptr<FILE, FileDeleter>(fopen(filename.c_str(), mode.c_str()));
    if (ret == nullptr) {
        std::ostringstream errOss;
        errOss << "Could not open file '" << filename << "'!";
        throw std::runtime_error(errOss.str());
    }
    return ret;
}

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_UTIL_H