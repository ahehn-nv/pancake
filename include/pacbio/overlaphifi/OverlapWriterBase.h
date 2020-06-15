// Author: Ivan Sovic

#ifndef PANCAKE_OVERLAPHIFI_OVERLAP_WRITER_BASE_H
#define PANCAKE_OVERLAPHIFI_OVERLAP_WRITER_BASE_H

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

class OverlapWriterBase
{
public:
    virtual ~OverlapWriterBase() {}

    virtual void Write(const OverlapPtr& ovl, const PacBio::Pancake::SeqDBReaderCached& targetSeqs,
                       const PacBio::Pancake::FastaSequenceId& querySeq, bool isFlipped) = 0;

    virtual void Write(const OverlapPtr& ovl,
                       const PacBio::Pancake::SeqDBReaderCachedBlock& targetSeqs,
                       const PacBio::Pancake::FastaSequenceCached& querySeq, bool isFlipped) = 0;

    static void PrintOverlapAsIPAOvl(FILE* fpOut, const OverlapPtr& ovl, const std::string& Aname,
                                     const std::string& Bname, bool writeIds, bool writeCigar);
    static void PrintOverlapAsM4(FILE* fpOut, const OverlapPtr& ovl, const std::string& Aname,
                                 const std::string& Bname, bool writeIds, bool writeCigar);
    static void PrintOverlapAsSAM(FILE* fpOut, const OverlapPtr& ovl, const char* seq,
                                  int64_t seqLen, const std::string& Aname,
                                  const std::string& Bname, bool writeIds, bool writeCigar);
    static std::string PrintOverlapAsM4(const OverlapPtr& ovl, const std::string& Aname,
                                        const std::string& Bname, bool writeIds, bool writeCigar);
};

constexpr char ConstexprTypeToChar(PacBio::BAM::CigarOperationType type)
{
    constexpr char lookup[11] = "MIDNSHP=XB";
    const int32_t x = static_cast<int>(type);
    return lookup[x];
}

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_OVERLAPHIFI_OVERLAP_WRITER_BASE_H
