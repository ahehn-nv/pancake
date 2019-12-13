// Authors: Ivan Sovic

#include "seqdb/SeqDBWorkflow.h"
#include <pbbam/FastaReader.h>
#include <seqdb/SeqDBWriterCompressed.h>
#include "seqdb/SeqDBSettings.h"

#include <string>
#include <vector>

#include <iostream>

namespace PacBio {
namespace Pancake {

int SeqDBWorkflow::Runner(const PacBio::CLI_v2::Results& options)
{
    SeqDBSettings settings{options};

    auto writer = PacBio::Pancake::CreateSeqDBWriterCompressed(
        settings.OutputPrefix, settings.BufferSize, settings.BlockSize);

    for (const auto& inFile : settings.InputFiles) {
        BAM::FastaReader inReader{inFile};
        BAM::FastaSequence record;
        while (inReader.GetNext(record)) {
            writer->AddSequence(record.Name(), record.Bases());
        }
    }

    return EXIT_SUCCESS;
}

}  // namespace Pancake
}  // namespace PacBio
