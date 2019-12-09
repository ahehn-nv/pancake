// Author: Ivan Sovic

#include <cstdlib>
#include <iostream>
#include <stdexcept>

#include <pbcopper/cli2/CLI.h>

#include <pacbio/Version.h>

#include "overlaphifi/OverlapHifiSettings.h"
#include "overlaphifi/OverlapHifiWorkflow.h"
#include "seqdb/SeqDBSettings.h"
#include "seqdb/SeqDBWorkflow.h"

PacBio::CLI_v2::MultiToolInterface CreateMultiInterface()
{
    PacBio::CLI_v2::MultiToolInterface mi{"pancake", "PacBio HiFi overlapper.",
                                          PacBio::Pancake::PancakeFormattedVersion()};

    // clang-format off
    mi.AddTools(
    {
        {"seqdb",
            PacBio::Pancake::SeqDBSettings::CreateCLI(),
           &PacBio::Pancake::SeqDBWorkflow::Runner},
        {"ovl-hifi",
            PacBio::Pancake::OverlapHifiSettings::CreateCLI(),
           &PacBio::Pancake::OverlapHifiWorkflow::Runner}
    });

    // clang-format on
    return mi;
}

int main(int argc, char* argv[])
{
    try {
        return PacBio::CLI_v2::Run(argc, argv, CreateMultiInterface());
    } catch (const std::exception& e) {
        std::cerr << "Pancake ERROR: " << e.what() << '\n';
        return EXIT_FAILURE;
    }
}