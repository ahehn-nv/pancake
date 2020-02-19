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
    OverlapWriter(const std::string& outFile, bool writeReverseOverlaps, bool writeIds);
    OverlapWriter(FILE* fpOut, bool writeReverseOverlaps, bool writeIds);
    ~OverlapWriter();

    void Write(const std::vector<OverlapPtr>& overlaps,
               const PacBio::Pancake::SeqDBReaderCached& targetSeqs,
               const PacBio::Pancake::FastaSequenceId& querySeq);
    void Write(const std::vector<OverlapPtr>& overlaps,
               const PacBio::Pancake::SeqDBReaderCachedBlock& targetSeqs,
               const PacBio::Pancake::FastaSequenceCached& querySeq);

    static void PrintOverlapAsM4(FILE* fpOut, const OverlapPtr& ovl, const std::string& Aname,
                                 const std::string& Bname, bool writeReverseOverlap, bool writeIds);
    static std::string PrintOverlapAsM4(const OverlapPtr& ovl, const std::string& Aname,
                                        const std::string& Bname, bool writeReverseOverlap,
                                        bool writeIds);

private:
    std::string outFile_;
    FILE* fpOut_ = NULL;
    bool shouldClose_ = false;
    bool writeReverseOverlaps_ = false;
    bool writeIds_ = false;
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_OVERLAPHIFI_OVERLAP_WRITER_H
