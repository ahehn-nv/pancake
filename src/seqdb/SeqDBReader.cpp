// Authors: Ivan Sovic

#include <pacbio/seqdb/CompressedSequence.h>
#include <pacbio/seqdb/SeqDBReader.h>
#include <pacbio/seqdb/Util.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>

namespace PacBio {
namespace Pancake {

SeqDBReader::SeqDBReader(std::shared_ptr<PacBio::Pancake::SeqDBIndexCache>& seqDBCache)
    : seqDBIndexCache_(seqDBCache)
{
    // Sanity check.
    if (seqDBIndexCache_->fileLines.empty())
        throw std::runtime_error("There are no file specifications in the input index file.");
    if (seqDBIndexCache_->seqLines.empty())
        throw std::runtime_error("There are no sequences in the input index file.");
}

SeqDBReader::~SeqDBReader() = default;

bool SeqDBReader::GetSequence(BAM::FastaSequence& record, int64_t seqId)
{
    // Clear the output.
    record.Bases("");
    record.Name("");

    // Sanity check for the sequence ID.
    if (seqId < 0 || seqId >= static_cast<int32_t>(seqDBIndexCache_->seqLines.size()))
        throw std::runtime_error("Invalid seqId out of bounds. seqId = " + std::to_string(seqId));

    // Find the sequence.
    auto it = seqDBIndexCache_->seqIdToOrdinalId.find(seqId);
    if (it == seqDBIndexCache_->seqIdToOrdinalId.end())
        throw std::runtime_error("Could not find seqId in the index. seqId = " +
                                 std::to_string(seqId));
    int32_t ordinalId = it->second;

    // Access the SequenceLine object.
    const auto& sl = seqDBIndexCache_->seqLines[ordinalId];

    // Load the data and decompress if required, and update fileHandler_->nextOrdinalId.
    LoadAndUnpackSequence_(record, fileHandler_, seqDBIndexCache_->fileLines,
                           seqDBIndexCache_->indexParentFolder, sl, ordinalId,
                           seqDBIndexCache_->compressionLevel);

    return true;
}

bool SeqDBReader::GetSequence(BAM::FastaSequence& record, const std::string& seqName)
{
    // Clear the output.
    record.Bases("");
    record.Name("");

    // Find the sequence.
    auto it = seqDBIndexCache_->headerToOrdinalId.find(seqName);
    if (it == seqDBIndexCache_->headerToOrdinalId.end()) return false;
    int32_t ordinalId = it->second;

    // Access the SequenceLine object.
    const auto& sl = seqDBIndexCache_->seqLines[ordinalId];

    // Load the data and decompress if required, and update fileHandler_->nextOrdinalId.
    LoadAndUnpackSequence_(record, fileHandler_, seqDBIndexCache_->fileLines,
                           seqDBIndexCache_->indexParentFolder, sl, ordinalId,
                           seqDBIndexCache_->compressionLevel);

    return true;
}

bool SeqDBReader::JumpTo(int64_t seqId)
{
    // Sanity check for the sequence ID.
    if (seqId < 0 || seqId >= static_cast<int32_t>(seqDBIndexCache_->seqLines.size()))
        throw std::runtime_error("Cannot JumpTo, invalid seqId out of bounds. seqId = " +
                                 std::to_string(seqId));

    // Find the sequence.
    auto it = seqDBIndexCache_->seqIdToOrdinalId.find(seqId);
    if (it == seqDBIndexCache_->seqIdToOrdinalId.end())
        throw std::runtime_error("Could not find seqId in the index. seqId = " +
                                 std::to_string(seqId));
    int32_t ordinalId = it->second;

    // Access the SequenceLine object.
    const auto& sl = seqDBIndexCache_->seqLines[ordinalId];

    // Jump to the correct file and offset, and update the fileHandler_->nextOrdinalId.
    AccessLocation_(fileHandler_, seqDBIndexCache_->fileLines, seqDBIndexCache_->indexParentFolder,
                    sl.fileId, ordinalId, sl.fileOffset);

    return true;
}

bool SeqDBReader::JumpTo(const std::string& seqName)
{
    // Find the sequence.
    auto it = seqDBIndexCache_->headerToOrdinalId.find(seqName);
    if (it == seqDBIndexCache_->headerToOrdinalId.end()) return false;
    int32_t ordinalId = it->second;

    // Access the SequenceLine object.
    const auto& sl = seqDBIndexCache_->seqLines[ordinalId];

    // Jump to the correct file and offset, and update the fileHandler_->nextOrdinalId.
    AccessLocation_(fileHandler_, seqDBIndexCache_->fileLines, seqDBIndexCache_->indexParentFolder,
                    sl.fileId, ordinalId, sl.fileOffset);

    return true;
}

bool SeqDBReader::GetNext(BAM::FastaSequence& record)
{
    // Clear the output.
    record.Bases("");
    record.Name("");

    // Sanity check for the sequence ID.
    if (fileHandler_.nextOrdinalId < 0) throw std::runtime_error("Invalid nextSeqId < 0.");

    // Can't go to the next sequence, we loaded all of them.
    if (fileHandler_.nextOrdinalId >= static_cast<int32_t>(seqDBIndexCache_->seqLines.size()))
        return false;

    // Access the SequenceLine object.
    const auto& sl = seqDBIndexCache_->seqLines[fileHandler_.nextOrdinalId];

    // Load the data and decompress if required, and update fileHandler_->nextOrdinalId.
    LoadAndUnpackSequence_(record, fileHandler_, seqDBIndexCache_->fileLines,
                           seqDBIndexCache_->indexParentFolder, sl, fileHandler_.nextOrdinalId,
                           seqDBIndexCache_->compressionLevel);

    return true;
}

bool SeqDBReader::GetNextBatch(std::vector<BAM::FastaSequence>& records, int64_t batchSize)
{
    records.clear();

    // Sanity check for the sequence ID.
    if (fileHandler_.nextOrdinalId < 0) throw std::runtime_error("Invalid nextSeqId < 0.");

    // Batch size < 0 loads everything as one batch.
    batchSize = (batchSize < 0) ? std::numeric_limits<int64_t>::max() : batchSize;

    int32_t numSeqLines = seqDBIndexCache_->seqLines.size();
    int64_t loadedBases = 0;
    while (fileHandler_.nextOrdinalId < numSeqLines) {
        // Access the SequenceLine object.
        const auto& sl = seqDBIndexCache_->seqLines[fileHandler_.nextOrdinalId];

        // Load the data and decompress if required.
        BAM::FastaSequence record;
        // Load the data and decompress if required, and update fileHandler_->nextOrdinalId.
        LoadAndUnpackSequence_(record, fileHandler_, seqDBIndexCache_->fileLines,
                               seqDBIndexCache_->indexParentFolder, sl, fileHandler_.nextOrdinalId,
                               seqDBIndexCache_->compressionLevel);

        // Store the record and increase the counters.
        records.emplace_back(record);
        loadedBases += static_cast<int64_t>(record.Bases().size());

        if (loadedBases >= batchSize) break;
    }

    if (records.empty()) return false;

    return true;
}

bool SeqDBReader::GetBlock(std::vector<BAM::FastaSequence>& records, int32_t blockId)
{
    records.clear();

    // Sanity check for the sequence ID.
    if (blockId < 0 || blockId >= static_cast<int32_t>(seqDBIndexCache_->blocks.size())) {
        std::ostringstream oss;
        oss << "Invalid blockId. blockId = " << blockId
            << ", blocks.size() = " << seqDBIndexCache_->blocks.size();
        throw std::runtime_error(oss.str());
    }

    const auto& block = seqDBIndexCache_->blocks[blockId];

    // Sanity check that the block's range is good.
    int32_t numSeqLines = seqDBIndexCache_->seqLines.size();
    if (block.start < 0 || block.end <= block.start || block.start >= numSeqLines ||
        block.end > numSeqLines) {
        std::ostringstream oss;
        oss << "Invalid block values: start = " << block.start << ", end = " << block.end
            << ", numSeqLines = " << numSeqLines << ".";
        throw std::runtime_error(oss.str());
    }

    fileHandler_.nextOrdinalId = block.start;

    while (fileHandler_.nextOrdinalId < block.end) {
        // Access the SequenceLine object.
        const auto& sl = seqDBIndexCache_->seqLines[fileHandler_.nextOrdinalId];

        // Load the data and decompress if required.
        BAM::FastaSequence record;
        // Load the data and decompress if required, and update fileHandler_->nextOrdinalId.
        LoadAndUnpackSequence_(record, fileHandler_, seqDBIndexCache_->fileLines,
                               seqDBIndexCache_->indexParentFolder, sl, fileHandler_.nextOrdinalId,
                               seqDBIndexCache_->compressionLevel);

        // Store the record.
        records.emplace_back(record);
    }

    if (records.empty()) return false;

    return true;
}

void SeqDBReader::AccessLocation_(OpenFileHandler& fileHandler,
                                  const std::vector<PacBio::Pancake::SeqDBFileLine>& fileLines,
                                  const std::string& indexParentFolder, int32_t fileId,
                                  int32_t nextOrdinalId, int64_t offset) const
{
    if (fileId < 0 || fileId >= static_cast<int32_t>(fileLines.size()))
        throw std::runtime_error("Invalid fileId value: " + std::to_string(fileId));
    if (offset < 0)
        throw std::runtime_error("Invalid file offset value: " + std::to_string(offset));

    const auto& fl = fileLines[fileId];

    if (fileHandler.fileId != fileId) {
        std::string actualPath = JoinPath(indexParentFolder, fl.filename);
        fileHandler.fp = PacBio::Pancake::OpenFile(actualPath.c_str(), "rb");
        fileHandler.fileId = fileId;
        fileHandler.pos = 0;
    }
    if (offset != fileHandler.pos) {
        int32_t rv = fseek(fileHandler.fp.get(), offset, SEEK_SET);
        if (rv) throw std::runtime_error("Could not fseek to position: " + std::to_string(offset));
        fileHandler.pos = offset;
    }
    fileHandler.nextOrdinalId = nextOrdinalId;
}

void SeqDBReader::LoadAndUnpackSequence_(
    BAM::FastaSequence& record, OpenFileHandler& fileHandler,
    const std::vector<PacBio::Pancake::SeqDBFileLine>& fileLines,
    const std::string& indexParentFolder, const SeqDBSequenceLine& sl, int32_t ordinalId,
    bool isCompressed) const
{

    // Jump to the correct file and offset.
    AccessLocation_(fileHandler, fileLines, indexParentFolder, sl.fileId, ordinalId + 1,
                    sl.fileOffset);

    if (isCompressed) {
        // Load and decompress the sequence.
        std::vector<uint8_t> data(sl.numBytes);
        int32_t n = fread(data.data(), sizeof(uint8_t), sl.numBytes, fileHandler.fp.get());
        if (n != sl.numBytes) {
            std::ostringstream oss;
            oss << "Could not read sequence of '" << sl.header << "'. Num bytes read: " << n
                << ", expected: " << sl.numBytes;
            throw std::runtime_error(oss.str());
        }
        std::string bases;
        DecompressSequence(data, sl.numBases, sl.ranges, bases);
        record.Name(sl.header);
        record.Bases(bases);

    } else {
        std::string bases(sl.numBytes, '\0');
        int32_t n = fread(&bases[0], sizeof(char), sl.numBytes, fileHandler.fp.get());
        if (n != sl.numBytes) {
            std::ostringstream oss;
            oss << "Could not read sequence of '" << sl.header << "'. Num bytes read: " << n
                << ", expected: " << sl.numBytes;
            throw std::runtime_error(oss.str());
        }
        record.Name(sl.header);
        record.Bases(bases);
    }
}

}  // namespace Pancake
}  // namespace PacBio
