// Author: Derek Barnett

#ifndef PANCAKE_WORKFLOW_H
#define PANCAKE_WORKFLOW_H

#include <pbcopper/cli2/CLI.h>

namespace PacBio {
namespace Pancake {

struct Workflow
{
    static int Runner(const PacBio::CLI_v2::Results& options);
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_WORKFLOW_H