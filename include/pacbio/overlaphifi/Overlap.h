// Author: Ivan Sovic

#ifndef PANCAKE_OVERLAPHIFI_OVERLAP_H
#define PANCAKE_OVERLAPHIFI_OVERLAP_H

#include <cstdint>
#include <memory>
#include <string>

namespace PacBio {
namespace Pancake {

enum class OverlapType
{
    Unknown,
    Internal,
    Contained,
    Contains,
    FivePrime,
    ThreePrime
};

class Overlap
{
public:
    int32_t Aid = -1;
    bool Arev = false;
    int32_t Astart = 0;
    int32_t Aend = 0;
    int32_t Alen = 0;

    int32_t Bid = -1;
    bool Brev = false;
    int32_t Bstart = 0;
    int32_t Bend = 0;
    int32_t Blen = 0;

    float Score = 0.0f;
    float Identity = 0.0f;
    int32_t EditDistance = -1;
    int32_t NumSeeds = -1;

    OverlapType Type = OverlapType::Unknown;

public:
    Overlap() = default;
    ~Overlap() = default;

    Overlap(int32_t _Aid, int32_t _Bid, float _Score, float _Identity, bool _Arev, int32_t _Astart,
            int32_t _Aend, int32_t _Alen, bool _Brev, int32_t _Bstart, int32_t _Bend, int32_t _Blen,
            int32_t _EditDistance, int32_t _NumSeeds, OverlapType _Type)
        : Aid(_Aid)
        , Arev(_Arev)
        , Astart(_Astart)
        , Aend(_Aend)
        , Alen(_Alen)
        , Bid(_Bid)
        , Brev(_Brev)
        , Bstart(_Bstart)
        , Bend(_Bend)
        , Blen(_Blen)
        , Score(_Score)
        , Identity(_Identity)
        , EditDistance(_EditDistance)
        , NumSeeds(_NumSeeds)
        , Type(_Type)
    {
    }

public:
    int32_t ASpan() const { return (Aend - Astart); }
    int32_t BSpan() const { return (Bend - Bstart); }
    int32_t AstartFwd() const { return (Arev ? (Alen - Aend) : Astart); }
    int32_t AendFwd() const { return (Arev ? (Alen - Astart) : Aend); }
    int32_t BstartFwd() const { return (Brev ? (Blen - Bend) : Bstart); }
    int32_t BendFwd() const { return (Brev ? (Blen - Bstart) : Bend); }

public:
    bool operator==(const Overlap& rhs) const
    {
        return Aid == rhs.Aid && Bid == rhs.Bid && Score == rhs.Score && Identity == rhs.Identity &&
               EditDistance == rhs.EditDistance && NumSeeds == rhs.NumSeeds && Type == rhs.Type &&
               Arev == rhs.Arev && Astart == rhs.Astart && Aend == rhs.Aend && Alen == rhs.Aend &&
               Brev == rhs.Brev && Bstart == rhs.Bstart && Bend == rhs.Bend && Blen == rhs.Bend;
    }
};

using OverlapPtr = std::unique_ptr<Overlap>;

inline std::unique_ptr<Overlap> createOverlap() { return std::unique_ptr<Overlap>(new Overlap()); }

inline std::unique_ptr<Overlap> createOverlap(int32_t Aid, int32_t Bid, float score, float identity,
                                              bool Arev, int32_t Astart, int32_t Aend, int32_t Alen,
                                              bool Brev, int32_t Bstart, int32_t Bend, int32_t Blen,
                                              int32_t EditDistance, int32_t NumSeeds,
                                              OverlapType Type)
{
    return std::unique_ptr<Overlap>(new Overlap(Aid, Bid, score, identity, Arev, Astart, Aend, Alen,
                                                Brev, Bstart, Bend, Blen, EditDistance, NumSeeds,
                                                Type));
}

inline std::unique_ptr<Overlap> createOverlap(const std::unique_ptr<Overlap>& ovl)
{
    return std::unique_ptr<Overlap>(new Overlap(
        ovl->Aid, ovl->Bid, ovl->Score, ovl->Identity, ovl->Arev, ovl->Astart, ovl->Aend, ovl->Alen,
        ovl->Brev, ovl->Bstart, ovl->Bend, ovl->Blen, ovl->EditDistance, ovl->NumSeeds, ovl->Type));
}

OverlapType DetermineOverlapType(const OverlapPtr& ovl, int32_t allowedDovetailDist);

void HeuristicExtendOverlapFlanks(OverlapPtr& ovl, int32_t allowedDist);

OverlapPtr HeuristicExtendOverlapFlanks(const OverlapPtr& ovl, int32_t allowedDist);

std::string OverlapTypeToString(const OverlapType& type);

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_OVERLAPHIFI_SEEDHIT_H