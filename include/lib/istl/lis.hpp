/*
 * lis.hpp
 *
 *      Author: Ivan Sovic
 *      GitHub: @isovic
 *      Copyright: Ivan Sovic, 2017, 2020
 *      Licence: MIT
 *
 * A generic implementation of the Longest Increasing Subsequence algorithm
 * which allows for custom data types, provided that a suitable comparison
 * function is given.
 */

#ifndef ISTL_LIS_H_
#define ISTL_LIS_H_

#include <cstdint>
#include <algorithm>
#include <vector>
#include <functional>

namespace istl {

template<class T>
std::vector<T> LIS(const std::vector<T> &points, int64_t begin, int64_t end,
                   std::function<bool(const T& a,
                                      const T& b)> compLessThan =
                                      [](const T& a, const T& b)
                                      { return a < b; } ) {
    /*
     * Based on the Python implementation here:
     * https://rosettacode.org/wiki/Longest_increasing_subsequence#Python
    */

    // Sanity check.
    if (points.size() == 0) {
        return {};
    }
    if (end < begin) {
        return {};
    }

    // Prepare the DP storage.
    const int64_t n = end - begin;
    std::vector<int64_t> dp(n + 1, 0);
    std::vector<int64_t> pred(n + 1, 0);
    int64_t len = 0;

    // Compute the LIS.
    for (int64_t i = 0; i < n; ++i) {
        int32_t low = 1;
        int32_t high = len;
        while (low <= high) {
            int32_t mid = (low + high) / 2;
            if (compLessThan(points[dp[mid] + begin], points[i + begin])) {
                low = mid + 1;
            } else {
                high = mid - 1;
            }
        }
        const int32_t newLen = low;
        pred[i] = dp[newLen - 1];
        dp[newLen] = i;
        if (newLen > len) {
            len = newLen;
        }
    }

    // Backtrack.
    std::vector<T> lis;
    lis.reserve(len);
    int64_t k = dp[len];
    for (int64_t i = (len - 1); i >= 0; --i) {
        lis.emplace_back(points[k + begin]);
        k = pred[k];
    }
    std::reverse(lis.begin(), lis.end());

    return lis;
}

}

#endif
