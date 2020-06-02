// Authors: Ivan Sovic

#include <pacbio/alignment/AlignmentTools.h>

#include <array>
#include <cstring>
#include <sstream>
#include <stdexcept>
#include <string>

#include <pbcopper/third-party/edlib.h>

namespace PacBio {
namespace Pancake {

PacBio::BAM::Cigar EdlibAlignmentToCigar(const unsigned char* aln, int32_t alnLen)
{
    if (alnLen <= 0) {
        return {};
    }

    // Edlib move codes: 0: '=', 1: 'I', 2: 'D', 3: 'X'
    std::array<PacBio::BAM::CigarOperationType, 4> opToCigar = {
        PacBio::BAM::CigarOperationType::SEQUENCE_MATCH, PacBio::BAM::CigarOperationType::INSERTION,
        PacBio::BAM::CigarOperationType::DELETION,
        PacBio::BAM::CigarOperationType::SEQUENCE_MISMATCH};

    PacBio::BAM::CigarOperationType prevOp = PacBio::BAM::CigarOperationType::UNKNOWN_OP;
    int32_t count = 0;
    PacBio::BAM::Cigar ret;
    for (int32_t i = 0; i <= alnLen; i++) {
        if (i == alnLen || (opToCigar[aln[i]] != prevOp &&
                            prevOp != PacBio::BAM::CigarOperationType::UNKNOWN_OP)) {
            ret.emplace_back(PacBio::BAM::CigarOperation(prevOp, count));
            count = 0;
        }
        if (i < alnLen) {
            prevOp = opToCigar[aln[i]];
            count += 1;
        }
    }
    return ret;
}

void EdlibAlignmentDiffCounts(const unsigned char* aln, int32_t alnLen, int32_t& numEq,
                              int32_t& numX, int32_t& numI, int32_t& numD)
{
    numEq = numX = numI = numD = 0;

    if (alnLen <= 0) {
        return;
    }

    for (int32_t i = 0; i < alnLen; i++) {
        switch (aln[i]) {
            case EDLIB_EDOP_MATCH:
                ++numEq;
                break;
            case EDLIB_EDOP_MISMATCH:
                ++numX;
                break;
            case EDLIB_EDOP_INSERT:
                ++numI;
                break;
            case EDLIB_EDOP_DELETE:
                ++numD;
                break;
            default:
                throw std::runtime_error("Unknown Edlib operation: " +
                                         std::to_string(static_cast<int32_t>(aln[i])));
                break;
        }
    }
}

void CigarDiffCounts(const PacBio::BAM::Cigar& cigar, int32_t& numEq, int32_t& numX, int32_t& numI,
                     int32_t& numD)
{
    numEq = numX = numI = numD = 0;

    if (cigar.empty()) {
        return;
    }

    for (const auto& op : cigar) {
        if (op.Type() == PacBio::BAM::CigarOperationType::SEQUENCE_MATCH) {
            numEq += op.Length();
        } else if (op.Type() == PacBio::BAM::CigarOperationType::SEQUENCE_MISMATCH) {
            numX += op.Length();
        } else if (op.Type() == PacBio::BAM::CigarOperationType::INSERTION) {
            numI += op.Length();
        } else if (op.Type() == PacBio::BAM::CigarOperationType::DELETION) {
            numD += op.Length();
        }
    }
}

void AppendToCigar(PacBio::BAM::Cigar& cigar, PacBio::BAM::CigarOperationType newOp, int32_t newLen)
{
    if (newLen == 0) {
        return;
    }
    if (cigar.empty() || newOp != cigar.back().Type()) {
        cigar.emplace_back(PacBio::BAM::CigarOperation(newOp, newLen));
    } else {
        cigar.back().Length(cigar.back().Length() + newLen);
    }
}

PacBio::BAM::Cigar ExpandMismatches(const char* query, int64_t queryLen, const char* target,
                                    int64_t targetLen, const PacBio::BAM::Cigar& cigar)
{
    PacBio::BAM::Cigar ret;
    if (cigar.size() <= 1) {
        return cigar;
    }
    int64_t queryPos = 0;
    int64_t targetPos = 0;
    int32_t lastAddedOp = -1;
    int32_t numCigarOps = cigar.size();

    for (int32_t i = 1; i < numCigarOps; ++i) {
        const auto& prevOp = cigar[i - 1];
        const auto& currOp = cigar[i];

        if (queryPos >= queryLen || targetPos >= targetLen) {
            std::ostringstream oss;
            oss << "Invalid CIGAR string: queryPos = " << queryPos << ", targetPos = " << targetPos
                << ", queryLen = " << queryLen << ", targetLen = " << targetLen;
            throw std::runtime_error(oss.str());
        }

        // Check if we have an INS+DEL or DEL+INS pair. If so, we'll convert them
        // into a single diagonal set of MATCH/MISMATCH operations, plus the left
        // hang and right hang indel operations.
        if ((prevOp.Type() == PacBio::BAM::CigarOperationType::INSERTION &&
             currOp.Type() == PacBio::BAM::CigarOperationType::DELETION) ||
            (prevOp.Type() == PacBio::BAM::CigarOperationType::DELETION &&
             currOp.Type() == PacBio::BAM::CigarOperationType::INSERTION)) {

            int32_t minLen = std::min(prevOp.Length(), currOp.Length());
            int32_t leftHang = static_cast<int32_t>(prevOp.Length()) - minLen;
            int32_t rightHang = static_cast<int32_t>(currOp.Length()) - minLen;

            AppendToCigar(ret, prevOp.Type(), leftHang);
            if (prevOp.Type() == PacBio::BAM::CigarOperationType::DELETION) {
                targetPos += leftHang;
            } else {
                queryPos += leftHang;
            }

            for (int32_t pos = 0; pos < minLen; ++pos) {
                if (query[queryPos] == target[targetPos]) {
                    AppendToCigar(ret, PacBio::BAM::CigarOperationType::SEQUENCE_MATCH, 1);
                } else {
                    AppendToCigar(ret, PacBio::BAM::CigarOperationType::SEQUENCE_MISMATCH, 1);
                }
                ++queryPos;
                ++targetPos;
            }

            AppendToCigar(ret, currOp.Type(), rightHang);
            if (currOp.Type() == PacBio::BAM::CigarOperationType::DELETION) {
                targetPos += rightHang;
            } else {
                queryPos += rightHang;
            }
            lastAddedOp = i;
            ++i;
            continue;
        } else {
            AppendToCigar(ret, prevOp.Type(), prevOp.Length());
            lastAddedOp = i - 1;
            if (prevOp.Type() != PacBio::BAM::CigarOperationType::DELETION) {
                queryPos += prevOp.Length();
            }
            if (prevOp.Type() != PacBio::BAM::CigarOperationType::INSERTION) {
                targetPos += prevOp.Length();
            }
        }
    }
    // Any remaining ops just get passed in.
    for (int32_t i = lastAddedOp + 1; i < numCigarOps; ++i) {
        AppendToCigar(ret, cigar[i].Type(), cigar[i].Length());
    }
    return ret;
}

void ValidateCigar(const char* query, int64_t queryLen, const char* target, int64_t targetLen,
                   const PacBio::BAM::Cigar& cigar, const std::string& label)
{
    int64_t queryPos = 0;
    int64_t targetPos = 0;
    int32_t numCigarOps = cigar.size();

    for (int32_t i = 0; i < numCigarOps; ++i) {
        const auto& op = cigar[i];

        if (queryPos > queryLen || targetPos > targetLen) {
            std::ostringstream oss;
            oss << "Invalid CIGAR string (global): "
                << "coordinates out of bounds! "
                << "queryPos = " << queryPos << ", targetPos = " << targetPos
                << ", offending CIGAR op: " << op.Length() << op.TypeToChar(op.Type())
                << ", queryLen = " << queryLen << ", targetLen = " << targetLen
                << ", CIGAR: " << cigar.ToStdString() << ", label: '" << label;
            throw std::runtime_error(oss.str());
        }

        if (op.Type() == PacBio::BAM::CigarOperationType::SEQUENCE_MATCH) {
            if ((queryPos + op.Length()) > queryLen || (targetPos + op.Length()) > targetLen) {
                std::ostringstream oss;
                oss << "Invalid CIGAR string (SEQUENCE_MATCH): "
                    << "coordinates out of bounds! "
                    << "queryPos = " << queryPos << ", targetPos = " << targetPos
                    << ", offending CIGAR op: " << op.Length() << op.TypeToChar(op.Type())
                    << ", queryLen = " << queryLen << ", targetLen = " << targetLen
                    << ", CIGAR: " << cigar.ToStdString() << ", label: '" << label;
                throw std::runtime_error(oss.str());
            }
            if (strncmp(query + queryPos, target + targetPos, op.Length()) != 0) {
                std::ostringstream oss;
                oss << "Invalid CIGAR string (SEQUENCE_MATCH): "
                    << "sequences are not equal even though they are delimited by a SEQUENCE_MATCH "
                       "operation! "
                    << "queryPos = " << queryPos << ", targetPos = " << targetPos
                    << ", offending CIGAR op: " << op.Length() << op.TypeToChar(op.Type())
                    << ", queryLen = " << queryLen << ", targetLen = " << targetLen
                    << ", CIGAR: " << cigar.ToStdString() << ", label: '" << label;
                throw std::runtime_error(oss.str());
            }
            queryPos += op.Length();
            targetPos += op.Length();
        } else if (op.Type() == PacBio::BAM::CigarOperationType::SEQUENCE_MISMATCH) {
            if ((queryPos + op.Length()) > queryLen || (targetPos + op.Length()) > targetLen) {
                std::ostringstream oss;
                oss << "Invalid CIGAR string (SEQUENCE_MISMATCH): "
                    << "coordinates out of bounds! "
                    << "queryPos = " << queryPos << ", targetPos = " << targetPos
                    << ", offending CIGAR op: " << op.Length() << op.TypeToChar(op.Type())
                    << ", queryLen = " << queryLen << ", targetLen = " << targetLen
                    << ", CIGAR: " << cigar.ToStdString() << ", label: '" << label;
                throw std::runtime_error(oss.str());
            }

            for (int64_t pos = 0; pos < static_cast<int64_t>(op.Length()); ++pos) {
                if (query[queryPos + pos] == target[targetPos + pos]) {
                    std::ostringstream oss;
                    oss << "Invalid CIGAR string (SEQUENCE_MISMATCH): "
                        << "sequences are equal even though they are delimited by a "
                           "SEQUENCE_MISMATCH operation! "
                        << "queryPos = " << queryPos << ", targetPos = " << targetPos
                        << ", offending CIGAR op: " << op.Length() << op.TypeToChar(op.Type())
                        << ", queryLen = " << queryLen << ", targetLen = " << targetLen
                        << ", CIGAR: " << cigar.ToStdString() << ", label: '" << label;
                    throw std::runtime_error(oss.str());
                }
            }
            queryPos += op.Length();
            targetPos += op.Length();
        } else if (op.Type() == PacBio::BAM::CigarOperationType::INSERTION ||
                   op.Type() == PacBio::BAM::CigarOperationType::SOFT_CLIP) {
            if ((queryPos + op.Length()) > queryLen) {
                std::ostringstream oss;
                oss << "Invalid CIGAR string (INSERTION): "
                    << "coordinates out of bounds! "
                    << "queryPos = " << queryPos << ", targetPos = " << targetPos
                    << ", offending CIGAR op: " << op.Length() << op.TypeToChar(op.Type())
                    << ", queryLen = " << queryLen << ", targetLen = " << targetLen
                    << ", CIGAR: " << cigar.ToStdString() << ", label: '" << label;
                throw std::runtime_error(oss.str());
            }
            queryPos += op.Length();
        } else if (op.Type() == PacBio::BAM::CigarOperationType::DELETION ||
                   op.Type() == PacBio::BAM::CigarOperationType::REFERENCE_SKIP) {
            if ((targetPos + op.Length()) > targetLen) {
                std::ostringstream oss;
                oss << "Invalid CIGAR string (DELETION): "
                    << "coordinates out of bounds! "
                    << "queryPos = " << queryPos << ", targetPos = " << targetPos
                    << ", offending CIGAR op: " << op.Length() << op.TypeToChar(op.Type())
                    << ", queryLen = " << queryLen << ", targetLen = " << targetLen
                    << ", CIGAR: " << cigar.ToStdString() << ", label: '" << label;
                throw std::runtime_error(oss.str());
            }
            targetPos += op.Length();
        } else if (op.Type() == PacBio::BAM::CigarOperationType::HARD_CLIP) {
            // Do nothing.
        } else {
            std::ostringstream oss;
            oss << "CIGAR operation '" << op.TypeToChar(op.Type())
                << "' not supported by the validation function.";
            throw std::runtime_error(oss.str());
        }
    }
}

void ExtractVariantString(const char* query, int64_t queryLen, const char* target,
                          int64_t targetLen, const PacBio::BAM::Cigar& cigar,
                          bool maskSimpleRepeats, bool maskHomopolymers,
                          std::string& retQueryVariants, std::string& retTargetVariants)
{
    int64_t queryPos = 0;
    int64_t targetPos = 0;
    int32_t numCigarOps = cigar.size();

    int32_t varStrQuerySize = 0;
    int32_t varStrTargetSize = 0;
    for (int32_t i = 0; i < numCigarOps; ++i) {
        const auto& op = cigar[i];
        if (op.Type() == PacBio::BAM::CigarOperationType::SEQUENCE_MISMATCH) {
            varStrQuerySize += op.Length();
            varStrTargetSize += op.Length();
        } else if (op.Type() == PacBio::BAM::CigarOperationType::INSERTION) {
            varStrQuerySize += op.Length();
        } else if (op.Type() == PacBio::BAM::CigarOperationType::DELETION) {
            varStrTargetSize += op.Length();
        }
    }

    std::string varStrQuery(varStrQuerySize, '0');
    std::string varStrTarget(varStrTargetSize, '0');
    int32_t varStrQueryPos = 0;
    int32_t varStrTargetPos = 0;

    for (int32_t i = 0; i < numCigarOps; ++i) {
        const auto& op = cigar[i];

        if (queryPos > queryLen || targetPos > targetLen) {
            std::ostringstream oss;
            oss << "Invalid CIGAR string (global): "
                << "coordinates out of bounds! "
                << "queryPos = " << queryPos << ", targetPos = " << targetPos
                << ", offending CIGAR op: " << op.Length() << op.TypeToChar(op.Type())
                << ", queryLen = " << queryLen << ", targetLen = " << targetLen
                << ", CIGAR: " << cigar.ToStdString();
            throw std::runtime_error(oss.str());
        }

        if (op.Type() == PacBio::BAM::CigarOperationType::SEQUENCE_MATCH) {
            // If it's a match, just move down the sequences.
            // Sanity check.
            if ((queryPos + op.Length()) > queryLen || (targetPos + op.Length()) > targetLen) {
                std::ostringstream oss;
                oss << "Invalid CIGAR string (SEQUENCE_MATCH): "
                    << "coordinates out of bounds! "
                    << "queryPos = " << queryPos << ", targetPos = " << targetPos
                    << ", offending CIGAR op: " << op.Length() << op.TypeToChar(op.Type())
                    << ", queryLen = " << queryLen << ", targetLen = " << targetLen
                    << ", CIGAR: " << cigar.ToStdString();
                throw std::runtime_error(oss.str());
            }
            // Move down.
            queryPos += op.Length();
            targetPos += op.Length();
        } else if (op.Type() == PacBio::BAM::CigarOperationType::SEQUENCE_MISMATCH) {
            // For a mismatch, include both alleles.
            // Sanity check.
            if ((queryPos + op.Length()) > queryLen || (targetPos + op.Length()) > targetLen) {
                std::ostringstream oss;
                oss << "Invalid CIGAR string (SEQUENCE_MISMATCH): "
                    << "coordinates out of bounds! "
                    << "queryPos = " << queryPos << ", targetPos = " << targetPos
                    << ", offending CIGAR op: " << op.Length() << op.TypeToChar(op.Type())
                    << ", queryLen = " << queryLen << ", targetLen = " << targetLen
                    << ", CIGAR: " << cigar.ToStdString();
                throw std::runtime_error(oss.str());
            }

            // Store the target allele first, then the query allele.
            for (int64_t pos = 0; pos < static_cast<int64_t>(op.Length()); ++pos) {
                varStrTarget[varStrTargetPos + pos] = target[targetPos + pos];
                varStrQuery[varStrQueryPos + pos] = query[queryPos + pos];
            }
            varStrQueryPos += op.Length();
            varStrTargetPos += op.Length();

            // Move down.
            queryPos += op.Length();
            targetPos += op.Length();
        } else if (op.Type() == PacBio::BAM::CigarOperationType::INSERTION) {
            // Sanity check.
            if ((queryPos + op.Length()) > queryLen) {
                std::ostringstream oss;
                oss << "Invalid CIGAR string (INSERTION): "
                    << "coordinates out of bounds! "
                    << "queryPos = " << queryPos << ", targetPos = " << targetPos
                    << ", offending CIGAR op: " << op.Length() << op.TypeToChar(op.Type())
                    << ", queryLen = " << queryLen << ", targetLen = " << targetLen
                    << ", CIGAR: " << cigar.ToStdString();
                throw std::runtime_error(oss.str());
            }

            // Add the query (insertion) bases.
            for (int64_t pos = 0; pos < static_cast<int64_t>(op.Length()); ++pos) {
                varStrQuery[varStrQueryPos + pos] = query[queryPos + pos];
            }
            varStrQueryPos += op.Length();

            // Move down.
            queryPos += op.Length();
        } else if (op.Type() == PacBio::BAM::CigarOperationType::DELETION) {
            // Sanity check.
            if ((targetPos + op.Length()) > targetLen) {
                std::ostringstream oss;
                oss << "Invalid CIGAR string (DELETION): "
                    << "coordinates out of bounds! "
                    << "queryPos = " << queryPos << ", targetPos = " << targetPos
                    << ", offending CIGAR op: " << op.Length() << op.TypeToChar(op.Type())
                    << ", queryLen = " << queryLen << ", targetLen = " << targetLen
                    << ", CIGAR: " << cigar.ToStdString();
                throw std::runtime_error(oss.str());
            }

            // Add the target (deletion) bases.
            for (int64_t pos = 0; pos < static_cast<int64_t>(op.Length()); ++pos) {
                varStrTarget[varStrTargetPos + pos] = target[targetPos + pos];
            }
            varStrTargetPos += op.Length();

            // Move down.
            targetPos += op.Length();

        } else if (op.Type() == PacBio::BAM::CigarOperationType::SOFT_CLIP) {
            // Sanity check.
            if ((queryPos + op.Length()) > queryLen) {
                std::ostringstream oss;
                oss << "Invalid CIGAR string (SOFT_CLIP): "
                    << "coordinates out of bounds! "
                    << "queryPos = " << queryPos << ", targetPos = " << targetPos
                    << ", offending CIGAR op: " << op.Length() << op.TypeToChar(op.Type())
                    << ", queryLen = " << queryLen << ", targetLen = " << targetLen
                    << ", CIGAR: " << cigar.ToStdString();
                throw std::runtime_error(oss.str());
            }

            // Move down.
            queryPos += op.Length();

        } else if (op.Type() == PacBio::BAM::CigarOperationType::REFERENCE_SKIP) {
            // Sanity check.
            if ((targetPos + op.Length()) > targetLen) {
                std::ostringstream oss;
                oss << "Invalid CIGAR string (REFERENCE_SKIP): "
                    << "coordinates out of bounds! "
                    << "queryPos = " << queryPos << ", targetPos = " << targetPos
                    << ", offending CIGAR op: " << op.Length() << op.TypeToChar(op.Type())
                    << ", queryLen = " << queryLen << ", targetLen = " << targetLen
                    << ", CIGAR: " << cigar.ToStdString();
                throw std::runtime_error(oss.str());
            }

            // Move down.
            targetPos += op.Length();

        } else {
            std::ostringstream oss;
            oss << "CIGAR operation '" << op.TypeToChar(op.Type())
                << "' not supported by ExtractVariantString function.";
            throw std::runtime_error(oss.str());
        }
    }

    std::swap(retQueryVariants, varStrQuery);
    std::swap(retTargetVariants, varStrTarget);
}

}  // namespace Pancake
}  // namespace PacBio
