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

#include <sstream>

namespace PacBio {
namespace Pancake {
namespace SeedDB {

static const int8_t MINIMIZER_FLAG_DEFAULT_FWD = 0x00;
static const int8_t MINIMIZER_FLAG_IS_REV = 0x01;
static const __int128 MINIMIZER_CODED_REV_BIT = (((__int128)1) << 32);

static const __int128 MINIMIZER_64bit_MASK = (((__int128)0x0FFFFFFFFFFFFFFFF));
static const __int128 MINIMIZER_32bit_MASK = (((__int128)0x0000000007FFFFFFF));
static const __int128 MINIMIZER_32bit_MASK_FULL = (((__int128)0x000000000FFFFFFFF));

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
    Seed(const __int128& codedKeypos)
        : key(DecodeKey(codedKeypos))
        , seqID(DecodeSeqId(codedKeypos))
        , pos(DecodePos(codedKeypos))
        , flag(DecodeIsRev(codedKeypos))
    {
    }

    bool IsRev() { return (flag & MINIMIZER_FLAG_IS_REV); }

    inline __int128 To128t()
    {
        __int128 ret = ((__int128)key) << 64;
        ret |= ((((__int128)(seqID)) << 1) & MINIMIZER_32bit_MASK) << 32;
        ret |= ((__int128)(pos)) & MINIMIZER_32bit_MASK;

        if (flag & MINIMIZER_FLAG_IS_REV) {
            ret |= MINIMIZER_CODED_REV_BIT;
        }

        return ret;
    }

    static inline __int128 Encode(uint64_t _key, int32_t _seqID, int32_t _pos, bool _isRev)
    {
        __int128 ret = ((__int128)_key) << 64;
        ret |= (((__int128)(_seqID & MINIMIZER_32bit_MASK)) << 33);
        ret |= ((__int128)(_pos & MINIMIZER_32bit_MASK));
        if (_isRev) {
            ret |= MINIMIZER_CODED_REV_BIT;
        }
        return ret;
    }

    static inline uint64_t DecodeKey(const __int128& seed)
    {
        return (uint64_t)((seed >> 64) & MINIMIZER_64bit_MASK);
    }

    static inline int32_t DecodePos(const __int128& seed)
    {
        return ((int32_t)(seed & MINIMIZER_32bit_MASK));
    }

    static inline int32_t DecodeSeqId(const __int128& seed)
    {
        return (int32_t)((seed >> 33) & MINIMIZER_32bit_MASK);
    }

    /*
     * Unlike DecodeSeqId, this method keeps the info about
     * of the strand still encoded in the return value.
     * The LSB is 0 for fwd, and 1 for rev.
     */
    static inline int32_t DecodeSeqIdWithRev(const __int128& seed)
    {
        return (int32_t)((seed >> 32) & MINIMIZER_32bit_MASK);
    }

    static inline bool DecodeIsRev(const __int128& seed)
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