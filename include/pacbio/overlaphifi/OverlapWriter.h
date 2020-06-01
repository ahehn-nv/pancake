// Author: Ivan Sovic

#ifndef PANCAKE_OVERLAPHIFI_OVERLAP_WRITER_H
#define PANCAKE_OVERLAPHIFI_OVERLAP_WRITER_H

#include <pacbio/overlaphifi/Overlap.h>
#include <pacbio/seqdb/FastaSequenceId.h>
#include <pacbio/seqdb/SeqDBReaderCached.h>
#include <pacbio/seqdb/SeqDBReaderCachedBlock.h>
#include <pacbio/seqdb/Util.h>
#include <cstdint>
#include <memory>
#include <unordered_map>
#include <vector>

namespace PacBio {
namespace Pancake {

class OverlapWriter
{
public:
    OverlapWriter(const std::string& outFile, bool writeReverseOverlaps,
                  int32_t allowedDovetailDist, bool writeIds, bool writeCigar);
    OverlapWriter(FILE* fpOut, bool writeReverseOverlaps, int32_t allowedDovetailDist,
                  bool writeIds, bool writeCigar);
    ~OverlapWriter();

    void Write(const std::vector<OverlapPtr>& overlaps,
               const PacBio::Pancake::SeqDBReaderCached& targetSeqs,
               const PacBio::Pancake::FastaSequenceId& querySeq);

    void Write(const std::vector<OverlapPtr>& overlaps,
               const PacBio::Pancake::SeqDBReaderCachedBlock& targetSeqs,
               const PacBio::Pancake::FastaSequenceCached& querySeq);

    static void PrintOverlapAsIPAOvl(FILE* fpOut, const OverlapPtr& ovl, const std::string& Aname,
                                     const std::string& Bname, bool writeIds, bool writeCigar);
    static void PrintOverlapAsM4(FILE* fpOut, const OverlapPtr& ovl, const std::string& Aname,
                                 const std::string& Bname, bool writeIds, bool writeCigar);
    static std::string PrintOverlapAsM4(const OverlapPtr& ovl, const std::string& Aname,
                                        const std::string& Bname, bool writeIds, bool writeCigar);

private:
    std::string outFile_;
    FILE* fpOut_ = NULL;
    bool shouldClose_ = false;
    bool writeReverseOverlaps_ = false;
    int32_t allowedDovetailDist_ = 0;
    bool writeIds_ = false;
    bool writeCigar_ = false;
};

constexpr char ConstexprTypeToChar(PacBio::BAM::CigarOperationType type)
{
    constexpr char lookup[11] = "MIDNSHP=XB";
    const int32_t x = static_cast<int>(type);
    return lookup[x];
}

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_OVERLAPHIFI_OVERLAP_WRITER_H
