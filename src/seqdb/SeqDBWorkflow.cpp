// Authors: Ivan Sovic

#include "seqdb/SeqDBWorkflow.h"
#include <pbbam/FastaReader.h>
#include <pbbam/FastqReader.h>
#include <seqdb/SeqDBWriter.h>
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

    auto isFasta = [](const std::string& fn) {

        return boost::algorithm::iends_with(fn, ".fasta") ||
               boost::algorithm::iends_with(fn, ".fasta.gz") ||
               boost::algorithm::iends_with(fn, ".fa") ||
               boost::algorithm::iends_with(fn, ".fa.gz");
    };
    auto isFastq = [](const std::string& fn) {
        return boost::algorithm::iends_with(fn, ".fastq") ||
               boost::algorithm::iends_with(fn, ".fastq.gz") ||
               boost::algorithm::iends_with(fn, ".fq") ||
               boost::algorithm::iends_with(fn, ".fq.gz");
    };

    for (const auto& inFile : settings.InputFiles) {
        if (isFasta(inFile)) {
            BAM::FastaReader inReader{inFile};
            BAM::FastaSequence record;
            while (inReader.GetNext(record)) {
                writer->AddSequence(record.Name(), record.Bases());
            }
        } else if (isFastq(inFile)) {
            BAM::FastqReader inReader{inFile};
            BAM::FastqSequence record;
            while (inReader.GetNext(record)) {
                writer->AddSequence(record.Name(), record.Bases());
            }
        } else {
            throw std::runtime_error("Unknown input file extension for file: '" + inFile + "'.");
        }
    }

    return EXIT_SUCCESS;
}

}  // namespace Pancake
}  // namespace PacBio
