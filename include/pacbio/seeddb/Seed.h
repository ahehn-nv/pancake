/*
 * seed.hpp
 *
 * This source file was obtained and adjusted from the Raptor
 * graph-based mapping tool codebase.
 *
 *  Created on: Oct 04, 2017
 *      Author: Ivan Sovic
 */

#ifndef SRC_MINIMIZER_INDEX2_MINIMIZER_H_
#define SRC_MINIMIZER_INDEX2_MINIMIZER_H_

#include <pacbio/util/CommonTypes.h>
#include <cstdint>
#include <sstream>

namespace PacBio {
namespace Pancake {
namespace SeedDB {

static const int8_t MINIMIZER_FLAG_DEFAULT_FWD = 0x00;
static const int8_t MINIMIZER_FLAG_IS_REV = 0x01;
static const PacBio::Pancake::Int128t MINIMIZER_CODED_REV_BIT =
    (((PacBio::Pancake::Int128t)1) << 32);

static const PacBio::Pancake::Int128t MINIMIZER_64bit_MASK =
    (((PacBio::Pancake::Int128t)0x0FFFFFFFFFFFFFFFF));
static const PacBio::Pancake::Int128t MINIMIZER_32bit_MASK =
    (((PacBio::Pancake::Int128t)0x0000000007FFFFFFF));
static const PacBio::Pancake::Int128t MINIMIZER_32bit_MASK_FULL =
    (((PacBio::Pancake::Int128t)0x000000000FFFFFFFF));

using SeedRaw = PacBio::Pancake::Int128t;

class Seed
{
public:
    Seed() : key(0), seqID(0), pos(0), flag(0) {}

    Seed(uint64_t _key, int32_t _seqID, int32_t _pos, bool _isRev)
        : key(_key)
        , seqID(_seqID)
        , pos(_pos)
        , flag((_isRev) ? MINIMIZER_FLAG_IS_REV : MINIMIZER_FLAG_DEFAULT_FWD)
    {
    }
    Seed(const PacBio::Pancake::Int128t& codedKeypos)
        : key(DecodeKey(codedKeypos))
        , seqID(DecodeSeqId(codedKeypos))
        , pos(DecodePos(codedKeypos))
        , flag(DecodeIsRev(codedKeypos))
    {
    }

    bool IsRev() { return (flag & MINIMIZER_FLAG_IS_REV); }

    inline PacBio::Pancake::Int128t To128t()
    {
        PacBio::Pancake::Int128t ret = ((PacBio::Pancake::Int128t)key) << 64;
        ret |= ((((PacBio::Pancake::Int128t)(seqID)) << 1) & MINIMIZER_32bit_MASK) << 32;
        ret |= ((PacBio::Pancake::Int128t)(pos)) & MINIMIZER_32bit_MASK;

        if (flag & MINIMIZER_FLAG_IS_REV) {
            ret |= MINIMIZER_CODED_REV_BIT;
        }

        return ret;
    }

    static inline PacBio::Pancake::Int128t Encode(uint64_t _key, int32_t _seqID, int32_t _pos,
                                                  bool _isRev)
    {
        PacBio::Pancake::Int128t ret = ((PacBio::Pancake::Int128t)_key) << 64;
        ret |= (((PacBio::Pancake::Int128t)(_seqID)) & MINIMIZER_32bit_MASK) << 33;
        ret |= ((PacBio::Pancake::Int128t)(_pos)) & MINIMIZER_32bit_MASK;
        if (_isRev) {
            ret |= MINIMIZER_CODED_REV_BIT;
        }
        return ret;
    }

    static inline uint64_t DecodeKey(const PacBio::Pancake::Int128t& seed)
    {
        return (uint64_t)((seed >> 64) & MINIMIZER_64bit_MASK);
    }

    static inline int32_t DecodePos(const PacBio::Pancake::Int128t& seed)
    {
        return ((int32_t)(seed & MINIMIZER_32bit_MASK));
    }

    static inline int32_t DecodeSeqId(const PacBio::Pancake::Int128t& seed)
    {
        return (int32_t)((seed >> 33) & MINIMIZER_32bit_MASK);
    }

    /*
     * Unlike DecodeSeqId, this method keeps the info about
     * of the strand still encoded in the return value.
     * The LSB is 0 for fwd, and 1 for rev.
     */
    static inline int32_t DecodeSeqIdWithRev(const PacBio::Pancake::Int128t& seed)
    {
        return (int32_t)((seed >> 32) & MINIMIZER_32bit_MASK);
    }

    static inline bool DecodeIsRev(const PacBio::Pancake::Int128t& seed)
    {
        return (seed & MINIMIZER_CODED_REV_BIT);
    }

    std::string Verbose() const
    {
        std::stringstream ss;
        ss << "pos = " << pos << ", seqID = " << seqID << ", flag = " << (int32_t)flag
           << ", key = " << key;
        return ss.str();
    }

    inline int32_t Compare(const Seed& other) const
    {
        if (key == other.key && seqID == other.seqID && pos == other.pos && flag == other.flag) {
            // Exact seed match.
            return 0;
        } else if (key != other.key) {
            // Keys are different, seeds do not match.
            return 1;
        }
        // Key is the same, the rest of the seed is different.
        return 2;
    }

    uint64_t key;
    int32_t seqID;
    int32_t pos;
    int8_t flag;
};

}  // namespace SeedDB
}  // namespace Pancake
}  // namespace PacBio

#endif
