// Author: Ivan Sovic

#ifndef PANCAKE_FILE_IO_H
#define PANCAKE_FILE_IO_H

#include <cstdint>
#include <string>
#include <unordered_set>
#include <vector>

namespace PacBio {
namespace Pancake {

enum class SequenceFormat
{
    Fasta,
    Fastq,
    Bam,
    Xml,
    Fofn,
    SeqDB,
    Unknown,
};

std::vector<std::string> LoadLinesToVector(const std::string& listPath);
std::unordered_set<std::string> LoadLinesToSet(const std::string& listPath);

bool FormatIsFasta(const std::string& fn);
bool FormatIsFastq(const std::string& fn);
bool FormatIsFofn(const std::string& fn);
bool FormatIsBam(const std::string& fn);
bool FormatIsXml(const std::string& fn);
bool FormatIsSeqDB(const std::string& fn);

SequenceFormat ParseFormat(const std::string& filename);

std::vector<std::pair<SequenceFormat, std::string>> ExpandInputFileList(
    const std::vector<std::string>& inFiles);

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_FILE_IO_H
