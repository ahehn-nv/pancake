// Authors:

#include "Workflow.h"

#include <pacbio/pancake/Settings.h>

namespace PacBio {
namespace Pancake {

int Workflow::Runner(const PacBio::CLI_v2::Results& options)
{
    Settings settings{options};

    //
    //
    // Ridiculously toy example, and not suggesting API/logic...
    // Just showing that component wiring & control flow go here in the 'runner'
    //
    //

    // OverlapDb overlapDb{settings};
    // OverlapFilters overlapFilters{settings};
    // if (shouldFilter)
    //     overlapDb.Filter(overlapFilters);

    // StringGraph graph{overlapDb, settings};
    // UnitigGraph unitigs{graph, settings};
    // ContigGraph contigs{unitigs, settings};

    // Reporter reporter{settings};
    // reporter.Print(contigs);

    return EXIT_SUCCESS;
}

}  // namespace Pancake
}  // namespace PacBio
