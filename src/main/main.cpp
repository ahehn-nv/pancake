// Author: Derek Barnett

#include <cstdlib>
#include <iostream>
#include <stdexcept>

#include <pbcopper/cli2/CLI.h>

#include <pacbio/pancake/Settings.h>
#include "Workflow.h"

int main(int argc, char* argv[])
{
    try {
        return PacBio::CLI_v2::Run(argc, argv, PacBio::Pancake::Settings::CreateCLI(),
                                   &PacBio::Pancake::Workflow::Runner);
    } catch (const std::runtime_error& e) {
        std::cerr << "Pancake ERROR: " << e.what() << '\n';
        return EXIT_FAILURE;
    }
}