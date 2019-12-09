// Authors: Ivan Sovic

#include "overlaphifi/OverlapHifiWorkflow.h"
#include "overlaphifi/OverlapHifiSettings.h"

namespace PacBio {
namespace Pancake {

int OverlapHifiWorkflow::Runner(const PacBio::CLI_v2::Results& options)
{
    OverlapHifiSettings settings{options};

    return EXIT_SUCCESS;
}

}  // namespace Pancake
}  // namespace PacBio
