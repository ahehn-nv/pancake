// Authors: Ivan Sovic

#include "seqdb/SeqDBWorkflow.h"
#include <pacbio/seqdb/SeqDBWriter.h>
#include <pacbio/util/FileIO.h>
#include <pbbam/BamReader.h>
#include <pbbam/DataSet.h>
#include <pbbam/FastaReader.h>
#include <pbbam/FastqReader.h>
#include <boost/algorithm/string/predicate.hpp>
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

    // Expand FOFNs and determine the formats of input files.
    std::vector<std::pair<SequenceFormat, std::string>> inputFiles =
        ExpandInputFileList(settings.InputFiles);

    for (const auto& inFilePair : inputFiles) {
        const auto& inFmt = inFilePair.first;
        const auto& inFile = inFilePair.second;
        if (inFmt == SequenceFormat::Fasta) {
            BAM::FastaReader inReader{inFile};
            BAM::FastaSequence record;
            while (inReader.GetNext(record)) {
                writer->AddSequence(record.Name(), record.Bases());
            }
        } else if (inFmt == SequenceFormat::Fastq) {
            BAM::FastqReader inReader{inFile};
            BAM::FastqSequence record;
            while (inReader.GetNext(record)) {
                writer->AddSequence(record.Name(), record.Bases());
            }
        } else if (inFmt == SequenceFormat::Bam) {
            BAM::BamReader inputBamReader{inFile};
            for (const auto& bam : inputBamReader)
                writer->AddSequence(bam.FullName(), bam.Sequence());
        } else {
            throw std::runtime_error("Unknown input file extension for file: '" + inFile + "'.");
        }
    }

    return EXIT_SUCCESS;
}

}  // namespace Pancake
}  // namespace PacBio
