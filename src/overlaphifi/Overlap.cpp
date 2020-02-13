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
