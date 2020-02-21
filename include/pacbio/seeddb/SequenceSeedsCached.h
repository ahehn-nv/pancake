// File Description
/// \file Sequence.h
/// \brief Defines the SequenceSeedsCached class. It contains all seeds for a given sequence,
//         tied with the sequence ID and length.
//
// Author: Ivan Sovic

#ifndef PANCAKE_SEQUENCE_SEEDS_CACHED_H
#define PANCAKE_SEQUENCE_SEEDS_CACHED_H

#include <cstdint>
#include <string>
#include <vector>

namespace PacBio {
namespace Pancake {

///
/// \brief The SequenceSeedsCached class hold a pointer for a given sequence, and
///        related data.
///
class SequenceSeedsCached
{
public:
    explicit SequenceSeedsCached(std::string name, const __int128* seeds, int64_t seedsSize,
                                 int64_t id);
    SequenceSeedsCached() = default;

    const std::string& Name() const;
    const __int128* Seeds() const;
    int64_t Size() const;
    int64_t Id() const;

    SequenceSeedsCached& Name(std::string name);
    SequenceSeedsCached& Seeds(const __int128* seeds);
    SequenceSeedsCached& Size(int64_t size);
    SequenceSeedsCached& Id(int64_t id);

private:
    std::string name_;
    const __int128* seeds_;
    int64_t size_;
    int64_t id_;
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_SEQUENCE_SEEDS_CACHED_H
