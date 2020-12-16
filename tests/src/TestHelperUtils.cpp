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

}  // namespace PancakeTests
}  // namespace PacBio
