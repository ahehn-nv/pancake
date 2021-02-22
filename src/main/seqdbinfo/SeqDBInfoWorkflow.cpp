// Authors: Ivan Sovic

#include "SeqDBInfoWorkflow.h"
#include "SeqDBInfoSettings.h"

#include <algorithm>
#include <cstdio>
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
    double totalLength = 0;
    int32_t numSeqs = 0;
    double lenMin = 0;
    double lenMax = 0;
    double lenAvg = 0.0;
    double lenMedian = 0.0;
    double nxAUC = 0.0;
    std::array<std::tuple<int32_t, double, int32_t>, 101> Nx;
    GenomicUnit unit;

    void Clear()
    {
        totalLength = numSeqs = lenMin = lenMax = 0;
        lenAvg = lenMedian = nxAUC = 0.0;
        unit = GenomicUnit::bp;
        for (size_t i = 0; i < Nx.size(); ++i) {
            Nx[i] = std::make_tuple(0, 0, 0);
        }
    }

    void ScaleLengthsByFactor(double factor)
    {
        totalLength *= factor;
        lenMin *= factor;
        lenMax *= factor;
        lenAvg *= factor;
        lenMedian *= factor;
        nxAUC *= factor;
        for (size_t i = 0; i < Nx.size(); ++i) {
            std::get<1>(Nx[i]) *= factor;
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

    // Compute the total length and Nx AUC.
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
    ret.nxAUC = (ret.totalLength == 0) ? 0.0 : ret.nxAUC / ret.totalLength;

    // Compute other basic stats.
    ret.unit = GenomicUnit::bp;
    ret.numSeqs = reverseSortedLengths.size();
    ret.lenMin = reverseSortedLengths.back();
    ret.lenMax = reverseSortedLengths.front();
    ret.lenAvg = (ret.numSeqs == 0) ? 0 : ret.totalLength / static_cast<double>(ret.numSeqs);
    ret.lenMedian = (ret.numSeqs == 0) ? 0
                                       : (reverseSortedLengths[ret.numSeqs / 2] +
                                          reverseSortedLengths[(ret.numSeqs - 1) / 2]) /
                                             2;

    // Compute the Nx graph for 1% increase in intervals.
    const int64_t nxGenomeSize = (genomeSize > 0) ? genomeSize : ret.totalLength;
    const double percStep = 0.01 * nxGenomeSize;
    double nextThreshold = 0;
    double subtotal = 0;
    int32_t x = 0;
    int32_t lastLen = reverseSortedLengths[0];
    int32_t lastLenId = 0;
    for (int32_t lenId = 0; lenId < static_cast<int32_t>(reverseSortedLengths.size()); ++lenId) {
        subtotal += reverseSortedLengths[lenId];
        while (subtotal >= nextThreshold && x < 101) {
            lastLen = reverseSortedLengths[lenId];
            lastLenId = lenId;
            ret.Nx[x] = std::make_tuple(x, lastLen, lastLenId + 1);
            nextThreshold += percStep;
            ++x;
        }
    }
    // Fill in any remaining bins.
    for (int32_t lenId = lastLenId; lenId < 101; ++lenId) {
        ret.Nx[x] = std::make_tuple(x, lastLen, lastLenId + 1);
    }
}

void ConvertStatsToUnit(SeqLengthStats& stats, const GenomicUnit& unit)
{
    const double factor1 = ConvertGenomicUnitToBpFactor(stats.unit);
    stats.ScaleLengthsByFactor(factor1);

    const double factor2 = ConvertBpToGenomicUnitFactor(unit);
    stats.ScaleLengthsByFactor(factor2);

    stats.unit = unit;
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
    result["unit"] = GenomicUnitToString(stats.unit);
    return result;
}

int SeqDBInfoWorkflow::Runner(const PacBio::CLI_v2::Results& options)
{
    SeqDBInfoSettings settings{options};

    // Load the SeqDB index cache.
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

    // Compute stats.
    SeqLengthStats stats;
    ComputeSeqLengthStats(lengths, 0, stats);

    // Convert to the desired unit.
    ConvertStatsToUnit(stats, settings.Unit);

    // Write the output.
    if (settings.HumanReadableOutput) {
        // clang-format off
        fprintf(stdout, "input\tunit\ttotal\tnum\tmin\tmax\tavg\tmedian\tAUC\tN10\tN10_n\tN25\tN25_n\tN50\tN50_n\tN75\tN75_n\tN90\tN90_n\tN100\tN100_n\n");

        fprintf(stdout, "%s\t%s\t%.2lf\t%d\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf"
                        "\t%.2lf\t%d\t%.2lf\t%d\t%.2lf\t%d\t%.2lf\t%d\t%.2lf\t%d\t%.2lf\t%d\n",
                        settings.InputSeqDB.c_str(),
                        GenomicUnitToString(stats.unit).c_str(),
                        stats.totalLength,
                        stats.numSeqs,
                        stats.lenMin,
                        stats.lenMax,
                        stats.lenAvg,
                        stats.lenMedian,
                        stats.nxAUC,
                        std::get<1>(stats.Nx[10]),
                        std::get<2>(stats.Nx[10]),
                        std::get<1>(stats.Nx[25]),
                        std::get<2>(stats.Nx[25]),
                        std::get<1>(stats.Nx[50]),
                        std::get<2>(stats.Nx[50]),
                        std::get<1>(stats.Nx[75]),
                        std::get<2>(stats.Nx[75]),
                        std::get<1>(stats.Nx[90]),
                        std::get<2>(stats.Nx[90]),
                        std::get<1>(stats.Nx[100]),
                        std::get<2>(stats.Nx[100])
        );
        // clang-format on

    } else {
        PacBio::JSON::Json statsJson = SeqLengthStatsToJson(stats);
        statsJson["total_bytes"] = totalBytes;
        // statsJson["input"] = settings.InputSeqDB;
        std::cout << statsJson.dump(4);
    }

    return EXIT_SUCCESS;
}

}  // namespace Pancake
}  // namespace PacBio
