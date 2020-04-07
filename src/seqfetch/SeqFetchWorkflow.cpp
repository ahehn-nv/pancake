// Authors: Ivan Sovic

#include "SeqFetchWorkflow.h"
#include "SeqFetchSettings.h"

#include <algorithm>
#include <iostream>
#include <string>
#include <tuple>
#include <unordered_set>
#include <vector>

#include <pacbio/seqdb/SeqDBIndexCache.h>
#include <pacbio/seqdb/SeqDBReader.h>
#include <pbbam/BamReader.h>
#include <pbbam/DataSet.h>
#include <pbbam/IndexedFastaReader.h>
#include <pbbam/IndexedFastqReader.h>
#include <pbbam/PbiFilterQuery.h>
#include <pbbam/PbiFilterTypes.h>
#include <boost/algorithm/string/predicate.hpp>

#include <pbcopper/logging/LogLevel.h>
#include <pbcopper/logging/Logging.h>

namespace PacBio {
namespace Pancake {

enum class SeqFetchInputFormat
{
    Fasta,
    Fastq,
    Bam,
    Xml,
    Fofn,
    SeqDB,
    Unknown,
};

std::vector<std::string> ParseVector(const std::string& listPath)
{
    std::vector<std::string> ret;
    std::ifstream ifs(listPath);
    if (ifs.is_open() == false) {
        throw std::runtime_error("Could not open file '" + listPath + "'!");
    }
    std::string line;
    while (std::getline(ifs, line)) {
        ret.emplace_back(line);
    }
    return ret;
}

std::unordered_set<std::string> ParseSet(const std::string& listPath)
{
    std::unordered_set<std::string> ret;
    std::ifstream ifs(listPath);
    if (ifs.is_open() == false) {
        throw std::runtime_error("Could not open file '" + listPath + "'!");
    }
    std::string line;
    while (std::getline(ifs, line)) {
        ret.emplace(line);
    }
    return ret;
}

bool FormatIsFasta(const std::string& fn)
{
    return boost::algorithm::iends_with(fn, ".fasta") ||
           boost::algorithm::iends_with(fn, ".fasta.gz") ||
           boost::algorithm::iends_with(fn, ".fa") || boost::algorithm::iends_with(fn, ".fa.gz");
}

bool FormatIsFastq(const std::string& fn)
{
    return boost::algorithm::iends_with(fn, ".fastq") ||
           boost::algorithm::iends_with(fn, ".fastq.gz") ||
           boost::algorithm::iends_with(fn, ".fq") || boost::algorithm::iends_with(fn, ".fq.gz");
}
bool FormatIsFofn(const std::string& fn) { return boost::algorithm::iends_with(fn, ".fofn"); }
bool FormatIsBam(const std::string& fn) { return boost::algorithm::iends_with(fn, ".bam"); }
bool FormatIsXml(const std::string& fn) { return boost::algorithm::iends_with(fn, ".xml"); }
bool FormatIsSeqDB(const std::string& fn) { return boost::algorithm::iends_with(fn, ".seqdb"); }

SeqFetchInputFormat ParseFormat(const std::string& filename)
{

    if (FormatIsFasta(filename)) {
        return SeqFetchInputFormat::Fasta;
    } else if (FormatIsFastq(filename)) {
        return SeqFetchInputFormat::Fastq;
    } else if (FormatIsBam(filename)) {
        return SeqFetchInputFormat::Bam;
    } else if (FormatIsXml(filename)) {
        return SeqFetchInputFormat::Xml;
    } else if (FormatIsFofn(filename)) {
        return SeqFetchInputFormat::Fofn;
    } else if (FormatIsSeqDB(filename)) {
        return SeqFetchInputFormat::SeqDB;
    }
    return SeqFetchInputFormat::Unknown;
}

std::vector<std::pair<SeqFetchInputFormat, std::string>> ExpandInputFileList(
    const std::vector<std::string>& inFiles)
{
    // Create the expanded list with loaded FOFNs, or BAM files from the XML.
    std::vector<std::pair<SeqFetchInputFormat, std::string>> retFiles;

    for (const auto& inFile : inFiles) {
        auto fmt = ParseFormat(inFile);

        if (fmt == SeqFetchInputFormat::Fofn) {
            std::vector<std::string> files = ParseVector(inFile);
            for (const auto& newFile : files) {
                auto newFmt = ParseFormat(newFile);
                retFiles.emplace_back(std::make_pair(newFmt, newFile));
            }

        } else if (fmt == SeqFetchInputFormat::Xml) {
            BAM::DataSet dataset{inFile};
            const auto& bamFiles = dataset.BamFiles();
            for (const auto& bamFile : bamFiles)
                retFiles.emplace_back(std::make_pair(SeqFetchInputFormat::Bam, bamFile.Filename()));

        } else {
            retFiles.emplace_back(std::make_pair(fmt, inFile));
        }
    }

    // Validate that all input formats are known.
    for (const auto& vals : retFiles) {
        const auto& fmt = std::get<0>(vals);
        const auto& file = std::get<1>(vals);
        if (fmt == SeqFetchInputFormat::Unknown) {
            throw std::runtime_error("Unknown extension of file: '" + file + "'.");
        }
    }

    return retFiles;
}

void WriteSeq(std::ostream& os, const std::string& seqName, const std::string& seq,
              const std::string& quals, char dummyQV, bool useProvidedQuals,
              const SeqFetchOutFormat& outFmt)
{
    if (outFmt == SeqFetchOutFormat::Fastq) {
        os << "@" << seqName << "\n" << seq << "\n+\n";
        if (useProvidedQuals) {
            os << quals << "\n";
        } else {
            os << std::string(seq.size(), dummyQV) << "\n";
        }
    } else if (outFmt == SeqFetchOutFormat::Fasta) {
        os << ">" << seqName << "\n" << seq << "\n";
    } else {
        throw std::runtime_error("Unknown output format!");
    }
}

void FetchFromFasta(std::shared_ptr<std::ostream>& osPtr, std::vector<std::string>& foundSeqs,
                    const std::string& inFile,
                    const std::unordered_set<std::string>& remainingToFind, const char dummyQV,
                    const PacBio::Pancake::SeqFetchOutFormat& outFormat)
{
    BAM::IndexedFastaReader reader{inFile};
    PBLOG_INFO << "Num sequences in file: " << reader.NumSequences();
    if (reader.NumSequences() == 0) {
        throw std::runtime_error("Index does not exist for file '" + inFile + "'!");
    }
    for (const auto& seqName : remainingToFind) {
        if (!reader.HasSequence(seqName)) continue;
        std::string seq = reader.Subsequence(seqName.c_str());
        std::string qual;
        foundSeqs.emplace_back(seqName);
        WriteSeq(*osPtr, seqName, seq, qual, dummyQV, false, outFormat);
    }
}

void FetchFromFastq(std::shared_ptr<std::ostream>& osPtr, std::vector<std::string>& foundSeqs,
                    const std::string& inFile,
                    const std::unordered_set<std::string>& remainingToFind, const char dummyQV,
                    const PacBio::Pancake::SeqFetchOutFormat& outFormat)
{
    BAM::IndexedFastqReader reader{inFile};
    PBLOG_INFO << "Num sequences in file: " << reader.NumSequences();
    if (reader.NumSequences() == 0) {
        throw std::runtime_error("Index does not exist for file '" + inFile + "'!");
    }
    for (const auto& seqName : remainingToFind) {
        if (!reader.HasSequence(seqName)) continue;
        int32_t seqLen = reader.SequenceLength(seqName);
        std::pair<std::string, PacBio::BAM::QualityValues> seqQualPair =
            reader.Subsequence(seqName, 0, seqLen);
        const auto& seq = std::get<0>(seqQualPair);
        std::string qual = std::get<1>(seqQualPair).Fastq();
        foundSeqs.emplace_back(seqName);
        WriteSeq(*osPtr, seqName, seq, qual, dummyQV, true, outFormat);
    }
}

void FetchFromBam(std::shared_ptr<std::ostream>& osPtr, std::vector<std::string>& foundSeqs,
                  const std::string& inFile, const std::unordered_set<std::string>& remainingToFind,
                  const char dummyQV, const PacBio::Pancake::SeqFetchOutFormat& outFormat)
{
    PacBio::BAM::BamFile bamFile(inFile);
    bamFile.EnsurePacBioIndexExists();

    std::vector<std::string> names;
    names.insert(names.end(), remainingToFind.begin(), remainingToFind.end());
    PacBio::BAM::PbiQueryNameFilter filter{names};
    PacBio::BAM::PbiFilterQuery query{filter, bamFile};
    for (auto record : query) {
        auto quals = record.Qualities().Fastq();
        WriteSeq(*osPtr, record.FullName(), record.Sequence(), quals, dummyQV, quals.size() > 0,
                 outFormat);
        foundSeqs.emplace_back(record.FullName());
    }
}

void FetchFromSeqDB(std::shared_ptr<std::ostream>& osPtr, std::vector<std::string>& foundSeqs,
                    const std::string& inFile,
                    const std::unordered_set<std::string>& remainingToFind, const char dummyQV,
                    const PacBio::Pancake::SeqFetchOutFormat& outFormat)
{
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBCache =
        PacBio::Pancake::LoadSeqDBIndexCache(inFile);

    PacBio::Pancake::SeqDBReader reader(seqDBCache);

    PacBio::Pancake::FastaSequenceId record;
    for (const auto& seqName : remainingToFind) {
        bool rv = reader.GetSequence(record, seqName);
        if (rv == false) {
            continue;
        }
        foundSeqs.emplace_back(seqName);
        WriteSeq(*osPtr, seqName, record.Bases(), std::string(), dummyQV, false, outFormat);
    }
}

int SeqFetchWorkflow::Runner(const PacBio::CLI_v2::Results& options)
{
    SeqFetchSettings settings{options};

    std::vector<std::pair<SeqFetchInputFormat, std::string>> inFiles =
        ExpandInputFileList(settings.InputFiles);
    std::unordered_set<std::string> seqNamesToFind = ParseSet(settings.InputFetchListFile);

    // If the input list is composed of sequence IDs, then the alias SeqDB
    // can be used to get the actual sequence names corresponding to those IDs.
    if (!settings.AliasSeqDBFile.empty()) {
        std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBCache =
            PacBio::Pancake::LoadSeqDBIndexCache(settings.AliasSeqDBFile);

        // The SeqDB is defined such that the ID matches the ordinal number of the seq line.
        std::vector<std::string> actualNames;
        for (const auto& idStr : seqNamesToFind) {
            int32_t id = std::stoi(idStr);
            const auto& sl = seqDBCache->GetSeqLine(id);
            actualNames.emplace_back(sl.header);
        }

        // Update the name list.
        seqNamesToFind.clear();
        seqNamesToFind.insert(actualNames.begin(), actualNames.end());
    }

    // Output file stream. If the file name is "-", then write to stdout.
    std::shared_ptr<std::ostream> osPtr(&std::cout, [](void*) {});
    if (settings.OutputFile.size() > 0 && settings.OutputFile != "-") {
        osPtr = std::shared_ptr<std::ostream>(new std::ofstream(settings.OutputFile));
        PBLOG_INFO << "Output is to file: " << settings.OutputFile;
    } else {
        PBLOG_INFO << "Output is to stdout.";
    }

    PBLOG_INFO << "Starting to fetch.";
    auto remainingToFind = seqNamesToFind;
    for (const auto& filePair : inFiles) {
        if (remainingToFind.empty()) {
            break;
        }

        const auto& inFmt = std::get<0>(filePair);
        const auto& inFile = std::get<1>(filePair);

        PBLOG_INFO << "Looking in file: " << inFile;

        std::vector<std::string> foundSeqs;

        if (inFmt == SeqFetchInputFormat::Fasta) {
            FetchFromFasta(osPtr, foundSeqs, inFile, remainingToFind, settings.DummyQV,
                           settings.OutputFormat);

        } else if (inFmt == SeqFetchInputFormat::Fastq) {
            FetchFromFastq(osPtr, foundSeqs, inFile, remainingToFind, settings.DummyQV,
                           settings.OutputFormat);

        } else if (inFmt == SeqFetchInputFormat::SeqDB) {
            FetchFromSeqDB(osPtr, foundSeqs, inFile, remainingToFind, settings.DummyQV,
                           settings.OutputFormat);

        } else if (inFmt == SeqFetchInputFormat::Bam) {
            FetchFromBam(osPtr, foundSeqs, inFile, remainingToFind, settings.DummyQV,
                         settings.OutputFormat);
        }

        // Remove the found sequences from the remaining set.
        for (const auto& found : foundSeqs) {
            remainingToFind.erase(found);
        }
    }

    PBLOG_INFO << "Done!";
    PBLOG_INFO << "Found sequences: " << (static_cast<int32_t>(seqNamesToFind.size()) -
                                          static_cast<int32_t>(remainingToFind.size()))
               << " / " << seqNamesToFind.size() << " .";

    if (settings.FailOnMissingQueries && remainingToFind.size() > 0) {
        std::ostringstream oss;
        oss << "Not all queries were found in the provided input files! Found sequences: "
            << (static_cast<int32_t>(seqNamesToFind.size()) -
                static_cast<int32_t>(remainingToFind.size()))
            << " / " << seqNamesToFind.size() << " .";
        throw std::runtime_error(oss.str());
    }

    return EXIT_SUCCESS;
}

}  // namespace Pancake
}  // namespace PacBio
