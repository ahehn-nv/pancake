// Authors: Ivan Sovic

#include "SeqDBInfoWorkflow.h"
#include "SeqDBInfoSettings.h"

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

#include <pacbio/pancake/SeqDBIndexCache.h>
#include <pbcopper/json/JSON.h>
#include <pbcopper/logging/LogLevel.h>
#include <pbcopper/logging/Logging.h>

namespace PacBio {
namespace Pancake {

struct SeqLengthStats
{
    int64_t totalLength = 0;
    int32_t numSeqs = 0;
    int32_t lenMin = 0;
    int32_t lenMax = 0;
    double lenAvg = 0.0;
    double lenMedian = 0.0;
    double nxAUC = 0.0;
    std::array<std::tuple<int32_t, int32_t, int32_t>, 101> Nx;

    void Clear()
    {
        totalLength = numSeqs = lenMin = lenMax = 0;
        lenAvg = lenMedian = nxAUC = 0.0;
        for (size_t i = 0; i < Nx.size(); ++i) {
            Nx[i] = std::make_tuple(0, 0, 0);
        }
    }
};

void ComputeSeqLengthStats(const std::vector<int32_t>& reverseSortedLengths, int64_t genomeSize,
                           SeqLengthStats& ret)
{
    ret.Clear();
    if (reverseSortedLengths.empty()) {
        return;
    }
    ret.totalLength = 0;
    ret.nxAUC = 0.0;
    int32_t prevVal = reverseSortedLengths.front();
    for (const int32_t& val : reverseSortedLengths) {
        if (val > prevVal) {
            throw std::runtime_error(
                "Input vector 'reverseSortedLengths' is not reverse sorted!\n");
        }
        ret.totalLength += val;
        ret.nxAUC += val * val;
    }
    ret.nxAUC = (ret.totalLength == 0) ? 0.0 : ret.nxAUC / static_cast<double>(ret.totalLength);
    ret.numSeqs = reverseSortedLengths.size();
    ret.lenMin = reverseSortedLengths.back();
    ret.lenMax = reverseSortedLengths.front();
    ret.lenAvg = (ret.numSeqs == 0)
                     ? 0
                     : static_cast<double>(ret.totalLength) / static_cast<double>(ret.numSeqs);
    ret.lenMedian = (ret.numSeqs == 0) ? 0
                                       : (reverseSortedLengths[ret.numSeqs / 2] +
                                          reverseSortedLengths[(ret.numSeqs - 1) / 2]) /
                                             2;

    const int64_t nxGenomeSize = (genomeSize > 0) ? genomeSize : ret.totalLength;
    const double percStep = 0.01 * nxGenomeSize;
    double nextThreshold = 0;  // percStep;
    double subtotal = 0;
    int32_t x = 0;
    for (int32_t lenId = 0; lenId < static_cast<int32_t>(reverseSortedLengths.size()); ++lenId) {
        subtotal += reverseSortedLengths[lenId];
        while (subtotal >= nextThreshold && x < 101) {
            ret.Nx[x] = std::make_tuple(x, reverseSortedLengths[lenId], lenId + 1);
            nextThreshold += percStep;
            ++x;
        }
    }
}

PacBio::JSON::Json SeqLengthStatsToJson(const SeqLengthStats& stats)
{
    auto result = PacBio::JSON::Json::object();
    result["total_length"] = stats.totalLength;
    result["num_seqs"] = stats.numSeqs;
    result["len_min"] = stats.lenMin;
    result["len_max"] = stats.lenMax;
    result["len_avg"] = stats.lenAvg;
    result["len_median"] = stats.lenMedian;
    result["AUC"] = stats.nxAUC;
    result["Nx"] = stats.Nx;
    return result;
}

int SeqDBInfoWorkflow::Runner(const PacBio::CLI_v2::Results& options)
{
    SeqDBInfoSettings settings{options};

    PBLOG_INFO << "Loading the SeqDB.";
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBCache =
        PacBio::Pancake::LoadSeqDBIndexCache(settings.InputSeqDB);

    // Collect and sort the length info.
    std::vector<int32_t> lengths(seqDBCache->seqLines.size(), 0);
    int64_t totalBytes = 0;
    for (size_t i = 0; i < seqDBCache->seqLines.size(); ++i) {
        const SeqDBSequenceLine& sl = seqDBCache->seqLines[i];
        lengths[i] = sl.numBases;
        totalBytes += sl.numBytes;
    }
    std::sort(lengths.begin(), lengths.end(), [](const auto& a, const auto& b) { return a > b; });

    SeqLengthStats stats;
    ComputeSeqLengthStats(lengths, 0, stats);

    PacBio::JSON::Json statsJson = SeqLengthStatsToJson(stats);
    statsJson["total_bytes"] = totalBytes;

    std::cout << statsJson.dump(4) << "\n";

    PBLOG_INFO << "Done!";

    return EXIT_SUCCESS;
}

}  // namespace Pancake
}  // namespace PacBio
