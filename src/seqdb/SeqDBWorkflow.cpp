// Authors: Ivan Sovic

#include "seqdb/SeqDBWorkflow.h"
#include <pbbam/FastaReader.h>
#include <seqdb/SeqDBWriter.h>
#include "seqdb/SeqDBSettings.h"

#include <string>
#include <vector>

#include <iostream>

namespace PacBio {
namespace Pancake {

int SeqDBWorkflow::Runner(const PacBio::CLI_v2::Results& options)
{
    SeqDBSettings settings{options};

    auto writer = PacBio::Pancake::CreateSeqDBWriter(settings.OutputPrefix,
                                                     settings.CompressionLevel, settings.BufferSize,
                                                     settings.BlockSize, settings.SplitBlocks);

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
