// Author: Ivan Sovic

#ifndef PANCAKE_TEST_HELPER_UTILS_H
#define PANCAKE_TEST_HELPER_UTILS_H

#include <pbbam/FastaSequence.h>

#include <string>
#include <vector>

namespace PacBio {
namespace PancakeTests {

std::vector<PacBio::BAM::FastaSequence> HelperLoadFasta(const std::string& inFasta);
std::vector<std::string> HelperLoadFastaAsStringVector(const std::string& inFasta);
std::string HelperLoadFastaAsString(const std::string& inFasta);

}  // namespace PancakeTests
}  // namespace PacBio

#endif  // PANCAKE_TEST_HELPER_UTILS_H