// Author: Ivan Sovic

#ifndef PANCAKE_RLE_H
#define PANCAKE_RLE_H

#include <cstdint>
#include <string>
#include <vector>

namespace PacBio {
namespace Pancake {

/// \brief Compresses the homopolymers from the input sequence. Returns (via
///         parameters) the compressed sequence, and the run lengths for each base.
void RunLengthEncoding(const std::string& seq, std::string& encodedSeq,
                       std::vector<int32_t>& runLengths);

/// \brief Compresses the homopolymers from the input sequence. Returns (via
///         parameters) the compressed sequence, and the run lengths for each base.
void RunLengthEncoding(const char* seq, int64_t seqLen, std::string& encodedSeq,
                       std::vector<int32_t>& runLengths);

/// \brief Compresses the homopolymers from the input sequence in-place.
int64_t RunLengthEncoding(char* seq, int64_t seqLen, std::vector<int32_t>& runLengths);

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_RLE_H