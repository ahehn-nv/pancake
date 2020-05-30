// Author: Ivan Sovic

#ifndef PANCAKE_ALIGNMENT_SES2_DISTANCE_BANDED_H
#define PANCAKE_ALIGNMENT_SES2_DISTANCE_BANDED_H

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <vector>

#include <pacbio/alignment/SesOptions.h>
#include <pacbio/alignment/SesResults.h>

// #define SES2_DEBUG

namespace PacBio {
namespace Pancake {
namespace Alignment {

template <SESAlignMode ALIGN_MODE, SESTrimmingMode TRIM_MODE>
SesResults SES2DistanceBanded(const char* query, size_t queryLen, const char* target,
                              size_t targetLen, int32_t maxDiffs, int32_t bandwidth)
{
    SesResults ret;

    if (queryLen == 0 || targetLen == 0) {
        ret.valid = true;
        return ret;
    }

    int32_t qlen = queryLen;
    int32_t tlen = targetLen;

    int32_t zero_offset = maxDiffs + 1;
    std::vector<int32_t> W(2 * maxDiffs + 3, MINUS_INF);    // Y for a diagonal k.
    std::vector<uint64_t> B;                                // Bitmask for trimming.
    std::vector<int32_t> M;                                 // Match count.

    if constexpr (TRIM_MODE == SESTrimmingMode::Enabled) {
        B.resize(2 * maxDiffs + 3, MINUS_INF);
        M.resize(2 * maxDiffs + 3, MINUS_INF);
    }

    std::vector<int32_t> u(2 * maxDiffs + 3, MINUS_INF);
    int32_t bandTolerance = bandwidth / 2 + 1;
    int32_t minK = 0;
    int32_t maxK = 0;
    int32_t best_u = 0;

    // Trimming related options.
    uint64_t b = 0;
    int32_t m = 0;
    const uint64_t C = 60;
    const uint64_t MASKC = static_cast<uint64_t>(1) << (C - 1);

    for (int32_t d = 0; d < maxDiffs; ++d) {
        ret.diffs = d;
        if ((maxK - minK) > bandwidth) {
            ret.valid = false;
            break;
        }

        int32_t ym = MINUS_INF;
        int32_t yc = MINUS_INF;
        int32_t yp = -1;

        W[zero_offset + minK - 1] = W[zero_offset + maxK + 1] = W[zero_offset + maxK + 0] = MINUS_INF;

        for (int32_t k = (minK); k <= (maxK); ++k) {
            int32_t kz = k + zero_offset;

            ym = yc;
            yc = yp;
            yp = W[kz + 1];

            int32_t minY = std::max(yc, std::max(ym, yp));

            int32_t y = MINUS_INF;

            if (ym == minY) {
                y = ym;
                if constexpr (TRIM_MODE == SESTrimmingMode::Enabled) {
                    m = M[kz - 1];
                    b = B[kz - 1];
                }
            } else if (yp == minY) {
                y = yp + 1; // Unlike 1986 paper, here we update y instead of x, so the +1 goes to the move to right (yp) instead of down (ym).
                if constexpr (TRIM_MODE == SESTrimmingMode::Enabled) {
                    m = M[kz + 1];
                    b = B[kz + 1];
                }
            } else {
                y = yc + ((yc < tlen) ? 1 : 0);
                if constexpr (TRIM_MODE == SESTrimmingMode::Enabled) {
                    m = M[kz];
                    b = B[kz];
                }
            }

            if constexpr (TRIM_MODE == SESTrimmingMode::Enabled) {
                if ((b & MASKC) != 0) {
                    --m;
                }
                b <<= 1;
            }

            int32_t x = y + k;

            while (x < qlen && y < tlen && query[x] == target[y]) {
                ++y;
                ++x;
                if constexpr (TRIM_MODE == SESTrimmingMode::Enabled) {
                    if ((b & MASKC) == 0) {
                        ++m;
                    }
                    b = (b << 1) | 1;
                }
            }

            W[kz] = y;
            if constexpr (TRIM_MODE == SESTrimmingMode::Enabled) {
                M[kz] = m;
                B[kz] = b;
            }

            u[kz] = y + k + y; // x + y = 2*y + k
            if (best_u <= u[kz]) {
                best_u = u[kz];
                ret.lastQueryPos = x;
                ret.lastTargetPos = y;
            }


            if constexpr(ALIGN_MODE == SESAlignMode::Global) {
                if (x >= qlen && y >= tlen) {
                    ret.valid = true;
                    ret.lastQueryPos = x;
                    ret.lastTargetPos = y;
                    break;
                }

            } else {
                if (x >= qlen || y >= tlen) {
                    ret.valid = true;
                    ret.lastQueryPos = x;
                    ret.lastTargetPos = y;
                    break;
                }
            }
        }

        if (ret.valid) {
            break;
        }

        int32_t newMinK = maxK;
        int32_t newMaxK = minK;
        for (int32_t k = (minK - 1); k <= (maxK + 1); ++k) {
            // Is there a bug here? Should this also have '&& u[k + zero_offset] <= (best_u + bandTolerance'?
            if (u[k + zero_offset] >= (best_u - bandTolerance)) {
                newMinK = std::min(k, newMinK);
                newMaxK = std::max(k, newMaxK);
            }
        }
        minK = newMinK - 1;
        maxK = newMaxK + 1;
    }

    return ret;
}
}
}
}

#endif  // PANCAKE_ALIGNMENT_SES_DISTANCE_BANDED_H
