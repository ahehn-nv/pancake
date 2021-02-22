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

inline double EnumToDouble(GenomicUnit unit)
{
    switch (unit) {
        case GenomicUnit::bp:
            return 1.0;
        case GenomicUnit::kbp:
            return 1'000.0;
        case GenomicUnit::Mbp:
            return 1'000'000.0;
        case GenomicUnit::Gbp:
            return 1'000'000'000.0;
        case GenomicUnit::Unknown:
        default:
            return 0.0;
    }
}

class GenomicUnitFromTo
{
public:  // until we update ScaleLengthByFactor() elsewhere
    double conversionFactor_;

public:
    GenomicUnitFromTo(GenomicUnit fromUnit, GenomicUnit toUnit)
        : conversionFactor_(EnumToDouble(fromUnit) / EnumToDouble(toUnit))
    {
    }
    double Convert(double fromValue) { return fromValue * conversionFactor_; }
    //int ConvertInt(int fromValue); // if needed
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_UTIL_GENOMIC_UNIT_H
