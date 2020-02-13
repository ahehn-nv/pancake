// Authors: Ivan Sovic

#include <pacbio/overlaphifi/OverlapWriter.h>

namespace PacBio {
namespace Pancake {

OverlapWriter::OverlapWriter(const std::string& outFile, bool writeReverseOverlaps, bool writeIds)
    : outFile_(outFile)
    , fpOut_(NULL)
    , shouldClose_(true)
    , writeReverseOverlaps_(writeReverseOverlaps)
    , writeIds_(writeIds)
{
    fpOut_ = fopen(outFile_.c_str(), "w");
}

OverlapWriter::OverlapWriter(FILE* fpOut, bool writeReverseOverlaps, bool writeIds)
    : outFile_("")
    , fpOut_(fpOut)
    , shouldClose_(false)
    , writeReverseOverlaps_(writeReverseOverlaps)
    , writeIds_(writeIds)
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
        PrintOverlapAsM4(fpOut_, ovl, querySeq.Name(), targetSeqs.GetSequence(ovl->Bid).Name(),
                         writeReverseOverlaps_, writeIds_);
    }
}

void OverlapWriter::PrintOverlapAsM4(FILE* fpOut, const OverlapPtr& ovl, const std::string& Aname,
                                     const std::string& Bname, bool writeReverseOverlap,
                                     bool writeIds)
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

    fprintf(fpOut, " %d %.2lf %d %d %d %d %d %d %d %d %s\n", static_cast<int32_t>(ovl->Score),
            100.0 * identity, static_cast<int32_t>(ovl->Arev), ovl->Astart, ovl->Aend, ovl->Alen,
            static_cast<int32_t>(tIsRev), tStart, tEnd, tLen, typeStr.c_str());

    if (writeReverseOverlap) {
        // The reverse overlap has the same coordinates as the normal orientation,
        // because all coordinates are in the FWD strand.
        // The Arev and Brev are intentionally kept in the normal orientation,
        // it's expected that the A-read is always fwd oriented.
        if (writeIds) {
            fprintf(fpOut, "%09d %09d", ovl->Bid, ovl->Aid);
        } else {
            fprintf(fpOut, "%s %s", Bname.c_str(), Aname.c_str());
        }
        OverlapType revType = (ovl->Type == OverlapType::FivePrime)
                                  ? OverlapType::ThreePrime
                                  : (ovl->Type == OverlapType::ThreePrime)
                                        ? OverlapType::FivePrime
                                        : (ovl->Type == OverlapType::Contained)
                                              ? OverlapType::Contains
                                              : (ovl->Type == OverlapType::Contains)
                                                    ? OverlapType::Contained
                                                    : ovl->Type;
        std::string revTypeStr = OverlapTypeToString(revType);
        fprintf(fpOut, " %d %.2lf %d %d %d %d %d %d %d %d %s\n", static_cast<int32_t>(ovl->Score),
                100.0 * identity, static_cast<int32_t>(ovl->Arev), tStart, tEnd, tLen,
                static_cast<int32_t>(tIsRev), ovl->Astart, ovl->Aend, ovl->Alen,
                revTypeStr.c_str());
    }
}

std::string OverlapWriter::PrintOverlapAsM4(const OverlapPtr& ovl, const std::string& Aname,
                                            const std::string& Bname, bool writeReverseOverlap,
                                            bool writeIds)
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

    oss << " " << static_cast<int32_t>(ovl->Score) << " " << idtBuff << " "
        << static_cast<int32_t>(ovl->Arev) << " " << ovl->Astart << " " << ovl->Aend << " "
        << ovl->Alen << " " << static_cast<int32_t>(tIsRev) << " " << tStart << " " << tEnd << " "
        << tLen << " " << typeStr.c_str();

    if (writeReverseOverlap) {
        // The reverse overlap has the same coordinates as the normal orientation,
        // because all coordinates are in the FWD strand.
        // The Arev and Brev are intentionally kept in the normal orientation,
        // it's expected that the A-read is always fwd oriented.
        if (writeIds) {
            oss << buffB << " " << buffA;
        } else {
            oss << Bname << " " << Aname;
        }
        OverlapType revType = (ovl->Type == OverlapType::FivePrime)
                                  ? OverlapType::ThreePrime
                                  : (ovl->Type == OverlapType::ThreePrime)
                                        ? OverlapType::FivePrime
                                        : (ovl->Type == OverlapType::Contained)
                                              ? OverlapType::Contains
                                              : (ovl->Type == OverlapType::Contains)
                                                    ? OverlapType::Contained
                                                    : ovl->Type;
        std::string revTypeStr = OverlapTypeToString(revType);
        oss << " " << static_cast<int32_t>(ovl->Score) << " " << 100.0 * identity << " "
            << static_cast<int32_t>(ovl->Arev) << " " << tStart << " " << tEnd << " " << tLen << " "
            << static_cast<int32_t>(tIsRev) << " " << ovl->Astart << " " << ovl->Aend << " "
            << ovl->Alen << " " << typeStr.c_str();
    }

    return oss.str();
}

}  // namespace Pancake
}  // namespace PacBio
