// Authors: Ivan Sovic

#include <pacbio/overlaphifi/OverlapWriter.h>

namespace PacBio {
namespace Pancake {

OverlapWriter::OverlapWriter(const std::string& outFile, bool writeReverseOverlaps,
                             int32_t allowedDovetailDist, bool writeIds, bool writeCigar)
    : outFile_(outFile)
    , fpOut_(NULL)
    , shouldClose_(true)
    , writeReverseOverlaps_(writeReverseOverlaps)
    , allowedDovetailDist_(allowedDovetailDist)
    , writeIds_(writeIds)
    , writeCigar_(writeCigar)
{
    fpOut_ = fopen(outFile_.c_str(), "w");
}

OverlapWriter::OverlapWriter(FILE* fpOut, bool writeReverseOverlaps, int32_t allowedDovetailDist,
                             bool writeIds, bool writeCigar)
    : outFile_("")
    , fpOut_(fpOut)
    , shouldClose_(false)
    , writeReverseOverlaps_(writeReverseOverlaps)
    , allowedDovetailDist_(allowedDovetailDist)
    , writeIds_(writeIds)
    , writeCigar_(writeCigar)
{
}

OverlapWriter::~OverlapWriter()
{
    if (shouldClose_) {
        fclose(fpOut_);
    }
}

void OverlapWriter::Write(const std::vector<OverlapPtr>& overlaps,
                          const PacBio::Pancake::SeqDBReaderCached& targetSeqs,
                          const PacBio::Pancake::FastaSequenceId& querySeq)
{
    for (const auto& ovl : overlaps) {
        // Don't look for the actual headers unless required. Saves the cost of a search.
        const auto& qName = writeIds_ ? "" : querySeq.Name();
        const auto& tName = writeIds_ ? "" : targetSeqs.GetSequence(ovl->Bid).Name();

        PrintOverlapAsM4(fpOut_, ovl, qName, tName, writeIds_, writeCigar_);

        if (writeReverseOverlaps_) {
            // Reverse overlaps cannot be collected in the Mapper class, because
            // the semantic of the Aread (which comes from the query DB) and the
            // Bread (which comes from the target DB) gets lost and it's impossible
            // to tell which DB Aread and Bread came from.
            auto newOvl = CreateFlippedOverlap(ovl, allowedDovetailDist_);
            PrintOverlapAsM4(fpOut_, newOvl, tName, qName, writeIds_, writeCigar_);
        }
    }
}

void OverlapWriter::Write(const std::vector<OverlapPtr>& overlaps,
                          const PacBio::Pancake::SeqDBReaderCachedBlock& targetSeqs,
                          const PacBio::Pancake::FastaSequenceCached& querySeq)
{
    for (const auto& ovl : overlaps) {
        // Don't look for the actual headers unless required. Saves the cost of a search.
        const auto& qName = writeIds_ ? "" : querySeq.Name();
        const auto& tName = writeIds_ ? "" : targetSeqs.GetSequence(ovl->Bid).Name();

        PrintOverlapAsM4(fpOut_, ovl, qName, tName, writeIds_, writeCigar_);

        if (writeReverseOverlaps_) {
            // Reverse overlaps cannot be collected in the Mapper class, because
            // the semantic of the Aread (which comes from the query DB) and the
            // Bread (which comes from the target DB) gets lost and it's impossible
            // to tell which DB Aread and Bread came from.
            auto newOvl = CreateFlippedOverlap(ovl, allowedDovetailDist_);
            PrintOverlapAsM4(fpOut_, newOvl, tName, qName, writeIds_, writeCigar_);
        }
    }
}

void OverlapWriter::PrintOverlapAsM4(FILE* fpOut, const OverlapPtr& ovl, const std::string& Aname,
                                     const std::string& Bname, bool writeIds, bool writeCigar)
{
    double identity = static_cast<double>(ovl->Identity);
    if (identity == 0.0 && ovl->EditDistance >= 0.0) {
        const double editDist = ovl->EditDistance;
        const double qSpan = ovl->ASpan();
        const double tSpan = ovl->BSpan();
        const double identityQ = (qSpan != 0) ? ((qSpan - editDist) / qSpan) : -2.0;
        const double identityT = (tSpan != 0) ? ((tSpan - editDist) / tSpan) : -2.0;
        identity = std::min(identityQ, identityT);
    }

    // The format specifies coordinates always in the FWD strand.
    int32_t tStart = ovl->BstartFwd();
    int32_t tEnd = ovl->BendFwd();
    const int32_t tIsRev = ovl->Brev;
    const int32_t tLen = ovl->Blen;
    std::string typeStr = OverlapTypeToString(ovl->Type);

    if (writeIds) {
        fprintf(fpOut, "%09d %09d", ovl->Aid, ovl->Bid);
    } else {
        fprintf(fpOut, "%s %s", Aname.c_str(), Bname.c_str());
    }

    fprintf(fpOut, " %d %.2lf %d %d %d %d %d %d %d %d %s", static_cast<int32_t>(ovl->Score),
            100.0 * identity, static_cast<int32_t>(ovl->Arev), ovl->Astart, ovl->Aend, ovl->Alen,
            static_cast<int32_t>(tIsRev), tStart, tEnd, tLen, typeStr.c_str());

    if (writeCigar) {
        if (ovl->Cigar.empty()) {
            fprintf(fpOut, " *");
        } else {
            fprintf(fpOut, " ");
            for (const auto& op : ovl->Cigar) {
                fprintf(fpOut, "%u%c", op.Length(), ConstexprTypeToChar(op.Type()));
            }
        }
    }

    fprintf(fpOut, "\n");
}

std::string OverlapWriter::PrintOverlapAsM4(const OverlapPtr& ovl, const std::string& Aname,
                                            const std::string& Bname, bool writeIds,
                                            bool writeCigar)
{
    double identity = static_cast<double>(ovl->Identity);
    if (identity == 0.0 && ovl->EditDistance >= 0.0) {
        const double editDist = ovl->EditDistance;
        const double qSpan = ovl->ASpan();
        const double tSpan = ovl->BSpan();
        const double identityQ = (qSpan != 0) ? ((qSpan - editDist) / qSpan) : -2.0;
        const double identityT = (tSpan != 0) ? ((tSpan - editDist) / tSpan) : -2.0;
        identity = std::min(identityQ, identityT);
    }

    // The format specifies coordinates always in the FWD strand.
    int32_t tStart = ovl->BstartFwd();
    int32_t tEnd = ovl->BendFwd();
    const int32_t tIsRev = ovl->Brev;
    const int32_t tLen = ovl->Blen;
    std::string typeStr = OverlapTypeToString(ovl->Type);

    std::ostringstream oss;
    char buffA[100], buffB[100];
    char idtBuff[100];
    sprintf(buffA, "%09d", ovl->Aid);
    sprintf(buffB, "%09d", ovl->Bid);
    sprintf(idtBuff, "%.2lf", 100.0 * identity);

    if (writeIds) {
        oss << buffA << " " << buffB;
    } else {
        oss << Aname << " " << Bname;
    }

    std::string cigar;
    std::string cigarSep;
    if (writeCigar) {
        cigarSep = " ";
        cigar = (ovl->Cigar.empty()) ? "*" : ovl->Cigar.ToStdString();
    }

    oss << " " << static_cast<int32_t>(ovl->Score) << " " << idtBuff << " "
        << static_cast<int32_t>(ovl->Arev) << " " << ovl->Astart << " " << ovl->Aend << " "
        << ovl->Alen << " " << static_cast<int32_t>(tIsRev) << " " << tStart << " " << tEnd << " "
        << tLen << " " << typeStr << cigarSep << cigar;

    return oss.str();
}

}  // namespace Pancake
}  // namespace PacBio
