// Authors: Ivan Sovic

#include <pacbio/seeddb/SeedDBIndexCache.h>
#include <pacbio/seqdb/Util.h>
#include <pbcopper/utility/StringUtils.h>
#include <limits>
#include <sstream>

namespace PacBio {
namespace Pancake {

std::unique_ptr<PacBio::Pancake::SeedDBIndexCache> LoadSeedDBIndexCache(
    const std::string& indexFilename)
{
    std::ifstream ifs(indexFilename);
    return LoadSeedDBIndexCache(ifs, indexFilename);
}

PacBio::Pancake::SeedDB::SeedDBParameters ParseSeedDBParams(const std::string& paramsStr)
{
    PacBio::Pancake::SeedDB::SeedDBParameters ret;

    auto params = PacBio::Utility::Split(paramsStr, ',');
    for (const auto& param : params) {
        if (param.empty()) continue;
        auto values = PacBio::Utility::Split(param, '=');
        if (values.size() != 2) {
            std::ostringstream oss;
            oss << "Parameter is not of form 'name=value'. Parameter: '" << param << "'.";
            throw std::runtime_error(oss.str());
        }
        if (values[0] == "k") {
            ret.KmerSize = std::stoi(values[1]);
        } else if (values[0] == "w") {
            ret.MinimizerWindow = std::stoi(values[1]);
        } else if (values[0] == "hpc") {
            ret.UseHPC = std::stoi(values[1]);
        } else if (values[0] == "hpc_len") {
            ret.MaxHPCLen = std::stoi(values[1]);
        } else if (values[0] == "rc") {
            ret.UseRC = std::stoi(values[1]);
        }
    }

    return ret;
}

std::unique_ptr<PacBio::Pancake::SeedDBIndexCache> LoadSeedDBIndexCache(
    std::istream& is, const std::string& indexFilename)
{
    auto cache = std::make_unique<PacBio::Pancake::SeedDBIndexCache>();

    cache->indexFilename = indexFilename;
    SplitPath(indexFilename, cache->indexParentFolder, cache->indexBasename);

    std::string line;
    char token;
    while (std::getline(is, line)) {
        if (line.empty()) continue;

        std::istringstream iss(line);
        iss >> token;

        SeedDBFileLine fl;
        SeedDBSeedsLine sl;
        SeedDBBlockLine bl;
        std::string paramsStr;
        int32_t ordinalId = 0;

        switch (token) {
            case 'V':
                iss >> cache->version;
                break;
            case 'P':
                iss >> paramsStr;
                cache->seedParams = ParseSeedDBParams(paramsStr);
                break;
            case 'F':
                iss >> fl.fileId >> fl.filename >> fl.numSequences >> fl.numBytes;
                cache->fileLines.emplace_back(fl);
                break;
            case 'S':
                iss >> sl.seqId >> sl.header >> sl.fileId >> sl.fileOffset >> sl.numBytes >>
                    sl.numBases >> sl.numSeeds;
                ordinalId = cache->seedLines.size();
                cache->seedLines.emplace_back(sl);
                // Add the new sequence to the lookups.
                cache->headerToOrdinalId[sl.header] = ordinalId;
                cache->seqIdToOrdinalId[sl.seqId] = ordinalId;
                break;
            case 'B':
                iss >> bl.blockId >> bl.startSeqId >> bl.endSeqId >> bl.numBytes;
                cache->blockLines.emplace_back(bl);
                break;
            default:
                std::ostringstream oss;
                oss << "Unknown token found when parsing the index: " << token;
                throw std::runtime_error(oss.str());
                break;
        }
    }

    if (cache->seedLines.empty())
        throw std::runtime_error("There are no sequences in the input index file: " +
                                 indexFilename);

    return cache;
}

const SeedDBSeedsLine& SeedDBIndexCache::GetSeedsLine(int32_t seqId) const
{
    // Find the seeds line.
    const auto it = seqIdToOrdinalId.find(seqId);
    if (it == seqIdToOrdinalId.end())
        throw std::runtime_error("Could not find seqId in the index. seqId = " +
                                 std::to_string(seqId));
    return seedLines[it->second];
}

const SeedDBSeedsLine& SeedDBIndexCache::GetSeedsLine(const std::string& seqName) const
{
    // Find the seeds line.
    const auto it = headerToOrdinalId.find(seqName);
    if (it == headerToOrdinalId.end())
        throw std::runtime_error(
            "Invalid sequence name, it does not exist in the SeedDB. seqName = " + seqName);
    return seedLines[it->second];
}

const SeedDBBlockLine& SeedDBIndexCache::GetBlockLine(int32_t blockId) const
{
    // Sanity check for the sequence ID.
    if (blockId < 0 || blockId >= static_cast<int32_t>(blockLines.size())) {
        std::ostringstream oss;
        oss << "Invalid blockId. blockId = " << blockId
            << ", blockLines.size() = " << blockLines.size();
        throw std::runtime_error(oss.str());
    }
    return blockLines[blockId];
}

const SeedDBFileLine& SeedDBIndexCache::GetFileLine(int32_t fileId) const
{
    // Sanity check.
    if (fileId < 0 || fileId >= static_cast<int32_t>(fileLines.size())) {
        std::ostringstream oss;
        oss << "Invalid fileId. fileId = " << fileId << ", fileLines.size() = " << fileLines.size();
        throw std::runtime_error(oss.str());
    }
    return fileLines[fileId];
}

std::ostream& operator<<(std::ostream& os, const PacBio::Pancake::SeedDBIndexCache& r)
{
    os << "V\t" << r.version << "\n";
    os << "P\tk=" << r.seedParams.KmerSize << ",w=" << r.seedParams.MinimizerWindow
       << ",hpc=" << r.seedParams.UseHPC << ",hpc_len=" << r.seedParams.MaxHPCLen
       << ",rc=" << r.seedParams.UseRC << "\n";
    for (const auto& fl : r.fileLines) {
        os << "F"
           << "\t" << fl.fileId << "\t" << fl.filename << "\t" << fl.numSequences << "\t"
           << fl.numBytes << "\n";
    }
    for (const auto& sl : r.seedLines) {
        os << "S"
           << "\t" << sl.seqId << "\t" << sl.header << "\t" << sl.fileId << "\t" << sl.fileOffset
           << "\t" << sl.numBytes << "\t" << sl.numBases << "\t" << sl.numSeeds << "\n";
    }
    for (const auto& bl : r.blockLines) {
        os << "B"
           << "\t" << bl.blockId << "\t" << bl.startSeqId << "\t" << bl.endSeqId << "\t"
           << bl.numBytes << "\n";
    }
    return os;
}

}  // namespace Pancake
}  // namespace PacBio
