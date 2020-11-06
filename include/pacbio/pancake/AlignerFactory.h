// Author: Ivan Sovic

#ifndef PANCAKE_ALIGNER_FACTORY_H
#define PANCAKE_ALIGNER_FACTORY_H

#include <pacbio/pancake/AlignerBase.h>
#include <pacbio/pancake/AlignerEdlib.h>
#include <pacbio/pancake/AlignerKSW2.h>
#include <pacbio/pancake/AlignerSES1.h>
#include <pacbio/pancake/AlignerSES2.h>
#include <pacbio/pancake/AlignerWFA.h>
#include <pacbio/pancake/AlignmentParameters.h>

namespace PacBio {
namespace Pancake {

enum class AlignerType
{
    KSW2,
    WFA,
    EDLIB,
    SES1,
    SES2,
};

std::string AlignerTypeToString(const AlignerType& alignerType);
AlignerType AlignerTypeFromString(const std::string& alignerType);
std::shared_ptr<AlignerBase> AlignerFactory(const AlignerType& alignerType,
                                            const AlignmentParameters& alnParams);

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_ALIGNER_FACTORY_H
