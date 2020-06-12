// Author: Ivan Sovic

#ifndef PANCAKE_OVERLAPHIFI_OVERLAP_WRITER_IPA_OVL_H
#define PANCAKE_OVERLAPHIFI_OVERLAP_WRITER_IPA_OVL_H

#include <pacbio/overlaphifi/Overlap.h>
#include <pacbio/overlaphifi/OverlapWriterBase.h>
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

class OverlapWriterIPAOvl : public OverlapWriterBase
{
public:
    OverlapWriterIPAOvl(FILE* fpOut, bool writeIds, bool writeCigar);
    ~OverlapWriterIPAOvl();

    void Write(const OverlapPtr& overlaps, const PacBio::Pancake::SeqDBReaderCached& targetSeqs,
               const PacBio::Pancake::FastaSequenceId& querySeq, bool isFlipped) override;

    void Write(const OverlapPtr& overlaps,
               const PacBio::Pancake::SeqDBReaderCachedBlock& targetSeqs,
               const PacBio::Pancake::FastaSequenceCached& querySeq, bool isFlipped) override;

private:
    std::string outFile_;
    FILE* fpOut_ = NULL;
    bool shouldClose_ = false;
    bool writeIds_ = false;
    bool writeCigar_ = false;
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_OVERLAPHIFI_OVERLAP_WRITER_IPA_OVL_H
