// Author: Ivan Sovic

#include "TestHelperUtils.h"
#include <pbbam/FastaReader.h>
#include <fstream>
#include <sstream>

namespace PacBio {
namespace PancakeTests {

std::vector<PacBio::BAM::FastaSequence> HelperLoadFasta(const std::string& inFasta)
{
    std::vector<PacBio::BAM::FastaSequence> ret;
    PacBio::BAM::FastaReader inReader{inFasta};
    PacBio::BAM::FastaSequence record;
    while (inReader.GetNext(record))
        ret.emplace_back(record);
    return ret;
}

std::string HelperLoadFastaAsString(const std::string& inFasta)
{
    std::ostringstream oss;
    auto records = HelperLoadFasta(inFasta);
    for (const auto& record : records)
        oss << ">" << record.Name() << "\n" << record.Bases() << "\n";
    return oss.str();
}

std::vector<std::string> HelperLoadFastaAsStringVector(const std::string& inFasta)
{
    std::vector<std::string> ret;
    PacBio::BAM::FastaReader inReader{inFasta};
    PacBio::BAM::FastaSequence record;
    while (inReader.GetNext(record)) {
        ret.emplace_back(record.Bases());
    }
    return ret;
}

std::vector<std::string> HelperLoadFile(const std::string& inFile)
{
    std::vector<std::string> ret;
    std::string line;
    std::ifstream ifs(inFile);
    if (ifs.is_open() == false) {
        throw std::runtime_error("Cannot open file " + inFile + " for reading!");
    }
    while (std::getline(ifs, line)) {
        ret.emplace_back(std::move(line));
    }
    return ret;
}

}  // namespace PancakeTests
}  // namespace PacBio
