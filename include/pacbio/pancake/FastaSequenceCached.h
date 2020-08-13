// Author: Ivan Sovic

#ifndef PANCAKE_FASTA_SEQUENCE_CACHED_H
#define PANCAKE_FASTA_SEQUENCE_CACHED_H

#include <cstdint>
#include <string>

namespace PacBio {
namespace Pancake {

class FastaSequenceCached
{
public:
    std::string name;
    const char* bases;
    int64_t size;
    int32_t id;

    const std::string& Name() const { return name; }
    const char* Bases() const { return bases; }
    int64_t Size() const { return size; }
    int32_t Id() const { return id; }
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_FASTA_SEQUENCE_CACHED_H