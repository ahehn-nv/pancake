// Authors: Ivan Sovic

#ifndef PANCAKE_SEQFETCH_SETTINGS_H
#define PANCAKE_SEQFETCH_SETTINGS_H

#include <pacbio/util/CommonTypes.h>
#include <cstdint>
#include <string>

#include <pbcopper/cli2/CLI.h>

namespace PacBio {
namespace Pancake {

enum class SeqFetchOutFormat
{
    Fasta,
    Fastq
};

struct SeqFetchSettings
{
    struct Defaults
    {
        static const char DummyQV = '!';
        static const bool FailOnMissingQueries = false;
    };

    std::string OutputFile;
    std::string InputFetchListFile;
    std::vector<std::string> InputFiles;
    SeqFetchOutFormat OutputFormat;
    char DummyQV = Defaults::DummyQV;
    std::string AliasSeqDBFile;
    bool FailOnMissingQueries = Defaults::FailOnMissingQueries;

    SeqFetchSettings();
    SeqFetchSettings(const PacBio::CLI_v2::Results& options);
    static PacBio::CLI_v2::Interface CreateCLI();
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_SEQFETCH_SETTINGS_H