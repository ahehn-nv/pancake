// Authors: Ivan Sovic

#include <pacbio/overlaphifi/Overlap.h>

namespace PacBio {
namespace Pancake {

OverlapType DetermineOverlapType(const OverlapPtr& ovl, int32_t allowedDovetailDist)
{
    int32_t leftHangA = ovl->AstartFwd();
    int32_t rightHangA = ovl->Alen - ovl->AendFwd();
    int32_t leftHangB = ovl->BstartFwd();
    int32_t rightHangB = ovl->Blen - ovl->BendFwd();

    if (ovl->Brev) {
        std::swap(leftHangB, rightHangB);
    }

    OverlapType ovlType = OverlapType::Internal;

    if (leftHangA <= allowedDovetailDist && rightHangA <= allowedDovetailDist) {
        /*
            This is a valid containment which would get picked up by the 5' rule
            if grace > 0. For this reason, the containment check must come first.
            Dovetail 5'.
                        left_a            right_a
                  >|/////////////|<          >|<
            A:     o-------------=============>
            B:                  o=============>
                               >|<         >||<
                               left_b     right_b

                               left_a     right_a
                               >|<         >||<
            A:                  o=============>
            B:     o-------------=============>
                  >|/////////////|<          >|<
                        left_b            right_b
        */
        ovlType = OverlapType::Contained;

    } else if (leftHangB <= allowedDovetailDist && rightHangB <= allowedDovetailDist) {
        ovlType = OverlapType::Contains;

    } else if (leftHangA <= allowedDovetailDist && rightHangB <= allowedDovetailDist) {
        /*
        Dovetail 5'.
                         left_a            right_a
                         >|//|<          >|///////|<
        A:                o--=============-------->
        B:     o-------------============->
              >|/////////////|<        >|/|<
                    left_b            right_b
        */
        ovlType = OverlapType::FivePrime;

    } else if (rightHangA <= allowedDovetailDist && leftHangB <= allowedDovetailDist) {
        /*
            Dovetail 5'.
                        left_a            right_a
                  >|/////////////|<        >|/|<
            A:     o-------------============->
            B:                o--=============-------->
                              >|//|<          >|///////|<
                              left_b            right_b
        */
        ovlType = OverlapType::ThreePrime;
    }

    return ovlType;
}

void HeuristicExtendOverlapFlanks(OverlapPtr& ovl, int32_t allowedDist)
{
    /// Note: The Overlap coordinates are internally represented in the strand
    /// of the overlap.

    int32_t leftHangA = ovl->Astart;
    int32_t rightHangA = ovl->Alen - ovl->Aend;
    int32_t leftHangB = ovl->Bstart;
    int32_t rightHangB = ovl->Blen - ovl->Bend;

    int32_t minLeft = std::min(leftHangA, leftHangB);
    int32_t minRight = std::min(rightHangA, rightHangB);

    int32_t leftFix = (minLeft > allowedDist) ? 0 : minLeft;
    int32_t rightFix = (minRight > allowedDist) ? 0 : minRight;

    ovl->Astart -= leftFix;
    ovl->Bstart -= leftFix;
    ovl->Aend += rightFix;
    ovl->Bend += rightFix;
}

OverlapPtr HeuristicExtendOverlapFlanks(const OverlapPtr& ovl, int32_t allowedDist)
{
    auto newOvl = createOverlap(ovl);
    HeuristicExtendOverlapFlanks(newOvl, allowedDist);
    return newOvl;
}

std::string OverlapTypeToString(const OverlapType& type)
{
    std::string ret = "x";
    switch (type) {
        case OverlapType::Unknown:
            ret = "x";
            break;
        case OverlapType::Contained:
            ret = "contained";
            break;
        case OverlapType::Contains:
            ret = "contains";
            break;
        case OverlapType::FivePrime:
            ret = "5";
            break;
        case OverlapType::ThreePrime:
            ret = "3";
            break;
        case OverlapType::Internal:
            ret = "u";
            break;
        default:
            ret = "x";
    }
    return ret;
}

}  // namespace Pancake
}  // namespace PacBio
