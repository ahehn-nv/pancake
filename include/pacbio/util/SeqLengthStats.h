// Author: Ivan Sovic

#ifndef PANCAKE_UTIL_SEQ_LENGTH_STATS_H
#define PANCAKE_UTIL_SEQ_LENGTH_STATS_H

#include <pacbio/util/GenomicUnit.h>
#include <pbcopper/json/JSON.h>
#include <memory>
#include <vector>

namespace PacBio {
namespace Pancake {

struct SeqLengthStats
{
    double totalLength = 0;
    int32_t numSeqs = 0;
    double lenMin = 0;
    double lenMax = 0;
    double lenAvg = 0.0;
    double lenMedian = 0.0;
    double nxAUC = 0.0;
    std::array<std::tuple<int32_t, double, int32_t>, 101> Nx;
    GenomicUnit unit;

    void Clear();
    void ChangeUnitOfLengths(GenomicUnit unit);

private:
    void ScaleLengthsByFactor_(double factor);
};

void ComputeSeqLengthStats(const std::vector<int32_t>& reverseSortedLengths, int64_t genomeSize,
                           SeqLengthStats& ret);

PacBio::JSON::Json SeqLengthStatsToJson(const SeqLengthStats& stats);

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_UTIL_SEQ_LENGTH_STATS_H