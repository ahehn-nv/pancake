// Author: Ivan Sovic

#ifndef PANCAKE_ALIGNMENT_TOOLS_H
#define PANCAKE_ALIGNMENT_TOOLS_H

#include <pacbio/alignment/DiffCounts.h>
#include <pbbam/Cigar.h>
#include <pbbam/CigarOperation.h>

namespace PacBio {
namespace Pancake {

PacBio::BAM::Cigar EdlibAlignmentToCigar(const unsigned char* aln, int32_t alnLen);

void EdlibAlignmentDiffCounts(const unsigned char* aln, int32_t alnLen, int32_t& numEq,
                              int32_t& numX, int32_t& numI, int32_t& numD);

void CigarDiffCounts(const PacBio::BAM::Cigar& cigar, int32_t& numEq, int32_t& numX, int32_t& numI,
                     int32_t& numD);

void AppendToCigar(PacBio::BAM::Cigar& cigar, PacBio::BAM::CigarOperationType newOp,
                   int32_t newLen);

PacBio::BAM::Cigar ExpandMismatches(const char* query, int64_t queryLen, const char* target,
                                    int64_t targetLen, const PacBio::BAM::Cigar& cigar);

void ValidateCigar(const char* query, int64_t queryLen, const char* target, int64_t targetLen,
                   const PacBio::BAM::Cigar& cigar, const std::string& label);

void ExtractVariantString(const char* query, int64_t queryLen, const char* target,
                          int64_t targetLen, const PacBio::BAM::Cigar& cigar, bool maskHomopolymers,
                          bool maskSimpleRepeats, std::string& retQueryVariants,
                          std::string& retTargetVariants, Alignment::DiffCounts& retDiffsPerBase,
                          Alignment::DiffCounts& retDiffsPerEvent);

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_ALIGNMENT_TOOLS_H