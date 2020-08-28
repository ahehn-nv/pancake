// Author: Ivan Sovic

#ifndef PANCAKE_OVERLAPHIFI_SECONDARY_H
#define PANCAKE_OVERLAPHIFI_SECONDARY_H

#include <lib/intervaltree/IntervalTree.h>
#include <pacbio/overlaphifi/Overlap.h>
#include <algorithm>
#include <cstdint>
#include <functional>
#include <memory>

using IntervalTreeInt32 =
    IntervalTree<int32_t, int32_t>;  // First: interval scalar type, Second: value type.
using IntervalVectorInt32 = IntervalTreeInt32::interval_vector;
using IntervalInt32 = IntervalTreeInt32::interval;

namespace PacBio {
namespace Pancake {

void FlagSecondaryAndSupplementary(std::vector<OverlapPtr>& overlaps, double allowedOverlapFraction,
                                   double minSecondaryScoreFraction);

void CreateRegionIntervalTrees(const std::vector<OverlapPtr>& overlaps,
                               std::function<bool(int32_t a)> CompPriority,
                               IntervalVectorInt32& queryIntervals, IntervalTreeInt32& queryTrees,
                               std::unordered_map<int32_t, IntervalVectorInt32>& targetIntervals,
                               std::unordered_map<int32_t, IntervalTreeInt32>& targetTrees);

bool CheckRegionSupplementary(const std::vector<OverlapPtr>& chainedRegions,
                              const OverlapPtr& currentRegion, IntervalTreeInt32& queryTrees,
                              std::unordered_map<int32_t, IntervalTreeInt32>& targetTrees,
                              double allowedOverlapFraction);

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_OVERLAPHIFI_SECONDARY_H
