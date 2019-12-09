// Authors: Ivan Sovic

#include "seqdb/SeqDBWorkflow.h"
#include "seqdb/SeqDBSettings.h"

namespace PacBio {
namespace Pancake {

int SeqDBWorkflow::Runner(const PacBio::CLI_v2::Results& options)
{
    SeqDBSettings settings{options};

    return EXIT_SUCCESS;
}

}  // namespace Pancake
}  // namespace PacBio
