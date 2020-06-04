// Authors: Ivan Sovic

#include <pacbio/overlaphifi/OverlapWriterBase.h>

namespace PacBio {
namespace Pancake {

void OverlapWriterBase::PrintOverlapAsIPAOvl(FILE* fpOut, const OverlapPtr& ovl,
                                             const std::string& Aname, const std::string& Bname,
                                             bool writeIds, bool writeCigar)
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

    // Format - 18 columns, space separated.
    //  Aid Bid score idt Arev Astart Aend Alen Brev Bstart Bend Blen ovl_type priority cigar Avars Bvars in_phase

    // The format specifies coordinates always in the FWD strand.
    int32_t tStart = ovl->BstartFwd();
    int32_t tEnd = ovl->BendFwd();
    const int32_t tIsRev = ovl->Brev;
    const int32_t tLen = ovl->Blen;
    std::string typeStr = OverlapTypeToString(ovl->Type);

    // First 12 columns are the same as in M4.
    if (writeIds) {
        fprintf(fpOut, "%09d %09d", ovl->Aid, ovl->Bid);
    } else {
        fprintf(fpOut, "%s %s", Aname.c_str(), Bname.c_str());
    }
    fprintf(fpOut, " %d %.6lf %d %d %d %d %d %d %d %d", static_cast<int32_t>(ovl->Score),
            100.0 * identity, static_cast<int32_t>(ovl->Arev), ovl->Astart, ovl->Aend, ovl->Alen,
            static_cast<int32_t>(tIsRev), tStart, tEnd, tLen);

    // Overlap type and priority. Priority is a placeholder for downstream tools.
    fprintf(fpOut, " %s *", typeStr.c_str());

    // Write the CIGAR only if specified, for speed.
    if (writeCigar) {
        if (ovl->Cigar.empty()) {
            fprintf(fpOut, " *");
        } else {
            fprintf(fpOut, " ");
            for (const auto& op : ovl->Cigar) {
                fprintf(fpOut, "%u%c", op.Length(), ConstexprTypeToChar(op.Type()));
            }
        }
    } else {
        fprintf(fpOut, " *");
    }

    // Write the A-read variant string, in the fwd orientation of the A-read.
    // Variant string is a list of variant bases for every non-match CIGAR operation.
    if (ovl->Avars.empty()) {
        fprintf(fpOut, " *");
    } else {
        if (ovl->Arev) {
            auto vars = Pancake::ReverseComplement(ovl->Avars, 0, ovl->Avars.size());
            fprintf(fpOut, " %s", vars.c_str());
        } else {
            fprintf(fpOut, " %s", ovl->Avars.c_str());
        }
    }

    // Write the B-read variant string, in the fwd orientation of the B-read.
    if (ovl->Bvars.empty()) {
        fprintf(fpOut, " *");
    } else {
        if (ovl->Brev) {
            auto vars = Pancake::ReverseComplement(ovl->Bvars, 0, ovl->Bvars.size());
            fprintf(fpOut, " %s", vars.c_str());
        } else {
            fprintf(fpOut, " %s", ovl->Bvars.c_str());
        }
    }

    // Write the "in-phase" info. At this point, everything is unphased.
    fprintf(fpOut, " u");
    fprintf(fpOut, "\n");
}

void OverlapWriterBase::PrintOverlapAsM4(FILE* fpOut, const OverlapPtr& ovl,
                                         const std::string& Aname, const std::string& Bname,
                                         bool writeIds, bool writeCigar)
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

void OverlapWriterBase::PrintOverlapAsSAM(FILE* fpOut, const OverlapPtr& ovl,
                                          const std::string& Aname, const std::string& Bname,
                                          bool writeIds, bool writeCigar)
{
    // double identity = static_cast<double>(ovl->Identity);
    // if (identity == 0.0 && ovl->EditDistance >= 0.0) {
    //     const double editDist = ovl->EditDistance;
    //     const double qSpan = ovl->ASpan();
    //     const double tSpan = ovl->BSpan();
    //     const double identityQ = (qSpan != 0) ? ((qSpan - editDist) / qSpan) : -2.0;
    //     const double identityT = (tSpan != 0) ? ((tSpan - editDist) / tSpan) : -2.0;
    //     identity = std::min(identityQ, identityT);
    // }

    // The format specifies coordinates always in the FWD strand.
    int32_t tStart = ovl->BstartFwd();
    int32_t tEnd = ovl->BendFwd();
    const int32_t tIsRev = ovl->Brev;
    const int32_t tLen = ovl->Blen;
    std::string typeStr = OverlapTypeToString(ovl->Type);
    int32_t flag = tIsRev ? 16 : 0;
    int32_t mapq = 60;
    std::string seq = "*";
    std::string qual = "*";

    // Query and target names, flag, pos and mapq.
    if (writeIds) {
        fprintf(fpOut, "%09d\t%d\t%09d\t%d\t%d", ovl->Aid, flag, ovl->Bid, tStart, mapq);
    } else {
        fprintf(fpOut, "%s\t%d\t%s\t%d\t%d", Aname.c_str(), flag, Bname.c_str(), tStart, mapq);
    }

    // Write the CIGAR only if specified, for speed.
    if (writeCigar) {
        if (ovl->Cigar.empty()) {
            fprintf(fpOut, "\t*");
        } else {
            fprintf(fpOut, "\t");
            for (const auto& op : ovl->Cigar) {
                fprintf(fpOut, "%u%c", op.Length(), ConstexprTypeToChar(op.Type()));
            }
        }
    } else {
        fprintf(fpOut, "\t*");
    }

    fprintf(fpOut, "\t*\t0\t0\t%s\t%s", seq.c_str(), qual.c_str());
    fprintf(fpOut, "\n");
}

std::string OverlapWriterBase::PrintOverlapAsM4(const OverlapPtr& ovl, const std::string& Aname,
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
