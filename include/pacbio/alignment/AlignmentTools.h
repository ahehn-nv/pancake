// Author: Ivan Sovic

#ifndef PANCAKE_ALIGNMENT_TOOLS_H
#define PANCAKE_ALIGNMENT_TOOLS_H

#include <pacbio/alignment/DiffCounts.h>
#include <pbbam/Cigar.h>
#include <pbbam/CigarOperation.h>

namespace PacBio {
namespace Pancake {

struct TrimmingInfo
{
    int32_t queryFront = 0;
    int32_t targetFront = 0;
    int32_t queryBack = 0;
    int32_t targetBack = 0;
};
inline bool operator==(const TrimmingInfo& lhs, const TrimmingInfo& rhs)
{
    return lhs.queryFront == rhs.queryFront && lhs.queryBack == rhs.queryBack &&
           lhs.targetFront == rhs.targetFront && lhs.targetBack == rhs.targetBack;
}

PacBio::BAM::Cigar EdlibAlignmentToCigar(const unsigned char* aln, int32_t alnLen);

void EdlibAlignmentDiffCounts(const unsigned char* aln, int32_t alnLen, int32_t& numEq,
                              int32_t& numX, int32_t& numI, int32_t& numD);

void CigarDiffCounts(const PacBio::BAM::Cigar& cigar, int32_t& numEq, int32_t& numX, int32_t& numI,
                     int32_t& numD);

Alignment::DiffCounts CigarDiffCounts(const PacBio::BAM::Cigar& cigar);

void AppendToCigar(PacBio::BAM::Cigar& cigar, PacBio::BAM::CigarOperationType newOp,
                   int32_t newLen);

PacBio::BAM::Cigar ExpandMismatches(const char* query, int64_t queryLen, const char* target,
                                    int64_t targetLen, const PacBio::BAM::Cigar& cigar);

void ValidateCigar(const char* query, int64_t queryLen, const char* target, int64_t targetLen,
                   const PacBio::BAM::Cigar& cigar, const std::string& label);

void ExtractVariantString(const char* query, int64_t queryLen, const char* target,
                          int64_t targetLen, const PacBio::BAM::Cigar& cigar, bool maskHomopolymers,
                          bool maskSimpleRepeats, bool maskHomopolymerSNPs,
                          bool maskHomopolymersArbitrary, std::string& retQueryVariants,
                          std::string& retTargetVariants, Alignment::DiffCounts& retDiffsPerBase,
                          Alignment::DiffCounts& retDiffsPerEvent);

Alignment::DiffCounts ComputeDiffCounts(const PacBio::BAM::Cigar& cigar,
                                        const std::string& queryVariants,
                                        const std::string& targetVariants,
                                        bool throwOnPartiallyMaskedIndels);

/// \brief For a given query position finds the corresponding target position based
///         on a provided CIGAR string.
///         Lineraly scans through all CIGAR operations to perform the mapping.
int32_t FindTargetPosFromCigar(const BAM::Cigar& cigar, int32_t queryPos);

void NormalizeAlignmentInPlace(std::string& queryAln, std::string& targetAln);

void ConvertCigarToM5(const char* query, int64_t queryLen, const char* target, int64_t targetLen,
                      const Data::Cigar& cigar, std::string& retQueryAln,
                      std::string& retTargetAln);

Data::Cigar ConvertM5ToCigar(const std::string& queryAln, const std::string& targetAln);

Data::Cigar NormalizeCigar(const char* query, int64_t queryLen, const char* target,
                           int64_t targetLen, const Data::Cigar& cigar);

bool TrimCigar(const PacBio::BAM::Cigar& cigar, int32_t windowSize, int32_t minMatches,
               bool clipOnFirstMatch, PacBio::BAM::Cigar& retTrimmedCigar,
               TrimmingInfo& retTrimming);

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_ALIGNMENT_TOOLS_H