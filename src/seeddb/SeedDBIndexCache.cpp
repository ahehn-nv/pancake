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
    FILE* fpIn = fopen(indexFilename.c_str(), "r");
    if (fpIn == NULL) {
        std::ostringstream oss;
        oss << "Could not open file '" << indexFilename << "' for reading!";
        throw std::runtime_error(oss.str());
    }
    auto result = LoadSeedDBIndexCache(fpIn, indexFilename);
    fclose(fpIn);
    return result;
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
    FILE* fpIn, const std::string& indexFilename)
{
    auto cache = std::make_unique<PacBio::Pancake::SeedDBIndexCache>();

    cache->indexFilename = indexFilename;
    SplitPath(indexFilename, cache->indexParentFolder, cache->indexBasename);

    char* line = NULL;
    size_t lineLen = 0;
    ssize_t numRead = 0;
    char buff[2000];  // Maximum string length (file names, headers).

    // Helper function to do cleanup of the C-based allocs.
    auto Cleanup = [&]() {
        if (line) {
            free(line);
            line = NULL;
            lineLen = 0;
        }
    };

    SeedDBFileLine fl;
    SeedDBSeedsLine sl;
    SeedDBBlockLine bl;
    int32_t ordinalId = 0;
    int32_t numReadItems = 0;
    size_t offset = 0;
    int32_t totalNumSeqs = 0;

    while ((numRead = getline(&line, &lineLen, fpIn)) != -1) {
        if (lineLen <= 0) {
            Cleanup();
            continue;
        }
        const char token = line[0];
        switch (token) {
            case 'V':
                numReadItems = sscanf(&line[1], "%s", buff);
                cache->version = buff;
                if (numReadItems != 1)
                    throw std::runtime_error("Problem parsing line: '" + std::string(line) + "'.");
                break;
            case 'P':
                // Find the first non-whitespace character.
                for (offset = 1; offset < lineLen && (offset == ' ' || offset == '\t'); ++offset) {
                }
                cache->seedParams = ParseSeedDBParams(std::string(line + offset));
                break;
            case 'F':
                numReadItems = sscanf(&line[1], "%d %s %d %lld", &(fl.fileId), buff,
                                      &(fl.numSequences), &(fl.numBytes));
                if (numReadItems != 4)
                    throw std::runtime_error("Problem parsing line: '" + std::string(line) + "'.");
                fl.filename = buff;
                cache->fileLines.emplace_back(fl);
                totalNumSeqs += fl.numSequences;
                cache->seedLines.reserve(totalNumSeqs);
                cache->headerToOrdinalId.reserve(totalNumSeqs);
                cache->seqIdToOrdinalId.reserve(totalNumSeqs);
                break;
            case 'S':
                numReadItems =
                    sscanf(&line[1], "%d %s %d %lld %lld %d %d", &(sl.seqId), buff, &(sl.fileId),
                           &(sl.fileOffset), &(sl.numBytes), &(sl.numBases), &(sl.numSeeds));
                if (numReadItems != 7)
                    throw std::runtime_error("Problem parsing line: '" + std::string(line) + "'.");
                sl.header = buff;
                ordinalId = cache->seedLines.size();
                cache->seedLines.emplace_back(sl);
                cache->headerToOrdinalId[sl.header] = ordinalId;
                cache->seqIdToOrdinalId[sl.seqId] = ordinalId;
                break;
            case 'B':
                numReadItems = sscanf(&line[1], "%d %d %d %lld", &(bl.blockId), &(bl.startSeqId),
                                      &(bl.endSeqId), &(bl.numBytes));
                if (numReadItems != 4)
                    throw std::runtime_error("Problem parsing line: '" + std::string(line) + "'.");
                cache->blockLines.emplace_back(bl);
                break;
            default:
                std::ostringstream oss;
                oss << "Unknown token found when parsing the index: " << token;
                throw std::runtime_error(oss.str());
                break;
        }
        Cleanup();
    }
    Cleanup();
    return cache;
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
