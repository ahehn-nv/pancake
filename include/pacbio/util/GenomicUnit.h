// Author: Ivan Sovic

#ifndef PANCAKE_UTIL_GENOMIC_UNIT_H
#define PANCAKE_UTIL_GENOMIC_UNIT_H

#include <pacbio/pancake/Lookups.h>
#include <cstdio>
#include <memory>
#include <sstream>
#include <vector>

namespace PacBio {
namespace Pancake {

enum class GenomicUnit
{
    bp,
    kbp,
    Mbp,
    Gbp,
    Unknown,
};

inline GenomicUnit GenomicUnitFromString(const std::string& val)
{
    if (val == "bp") {
        return GenomicUnit::bp;
    } else if (val == "kbp") {
        return GenomicUnit::kbp;
    } else if (val == "Mbp") {
        return GenomicUnit::Mbp;
    } else if (val == "Gbp") {
        return GenomicUnit::Gbp;
    }
    return GenomicUnit::Unknown;
}

inline std::string GenomicUnitToString(const GenomicUnit& unit)
{
    if (unit == GenomicUnit::bp) {
        return "bp";
    } else if (unit == GenomicUnit::kbp) {
        return "kbp";
    } else if (unit == GenomicUnit::Mbp) {
        return "Mbp";
    } else if (unit == GenomicUnit::Gbp) {
        return "Gbp";
    }
    return "Unknown";
}

inline double ConvertGenomicUnitToBpFactor(const GenomicUnit& sourceUnit)
{
    if (sourceUnit == GenomicUnit::bp) {
        return 1.0;
    } else if (sourceUnit == GenomicUnit::kbp) {
        return 1000.0;
    } else if (sourceUnit == GenomicUnit::Mbp) {
        return 1000000.0;
    } else if (sourceUnit == GenomicUnit::Gbp) {
        return 1000000000.0;
    }
    return 0.0;
}

inline double ConvertBpToGenomicUnitFactor(const GenomicUnit& targetUnit)
{
    if (targetUnit == GenomicUnit::bp) {
        return 1.0;
    } else if (targetUnit == GenomicUnit::kbp) {
        return 1.0 / 1000.0;
    } else if (targetUnit == GenomicUnit::Mbp) {
        return 1.0 / 1000000.0;
    } else if (targetUnit == GenomicUnit::Gbp) {
        return 1.0 / 1000000000.0;
    }
    return 0.0;
}

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_UTIL_GENOMIC_UNIT_H