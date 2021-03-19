/*
 * Minimizers.cpp
 *
 * This source file was obtained and adjusted from the Raptor
 * graph-based mapping tool codebase.
 *
 *  Created on: Oct 04, 2017
 *      Author: Ivan Sovic
 */

#include <pacbio/pancake/Lookups.h>
#include <pacbio/pancake/Minimizers.h>
#include <pacbio/pancake/Seed.h>
#include <pacbio/util/CommonTypes.h>
#include <deque>
#include <iostream>

namespace PacBio {
namespace Pancake {
namespace SeedDB {

using minkey_t = uint64_t;

const int32_t MAX_SPACING_IN_SEED = 32;
const int32_t MAX_WINDOW_BUFFER_SIZE = 512;

// This is a temporary class, just used for implementing the minimizer window.
class SpacedBuffer
{
public:
    SpacedBuffer()
    {
        for (size_t i = 0; i < MAX_WINDOW_BUFFER_SIZE; ++i) {
            winPosSet[i] = false;
        }

        for (size_t i = 0; i < seqBuffer.size(); ++i) {
            seqBuffer[i] = 0x0;
            seqBufferRC[i] = 0x0;
            numBasesIn[i] = 0;
        }
    }

    ~SpacedBuffer() = default;

    std::array<uint64_t, MAX_SPACING_IN_SEED>
        seqBuffer;  // Holds the current 2-bit seed representation.
    std::array<uint64_t, MAX_SPACING_IN_SEED>
        seqBufferRC;  // Holds the reverse complement 2-bit seed at the same position.
    std::array<int32_t, MAX_SPACING_IN_SEED> numBasesIn;  // Number of bases added to the buffer.

    Seed winBuff[MAX_WINDOW_BUFFER_SIZE];  // Define the new circular buffer for the window.
    bool winPosSet[MAX_WINDOW_BUFFER_SIZE];
    int32_t winBuffPos = 0;
    int32_t winBuffMinPos = -1;
    int32_t kmerSpan = 0;
    std::deque<int32_t> hpEvents;

    void Clear(int32_t winSize)
    {
        // buffer = 0;
        // buffer_rc = 0;
        // numBasesIn = 0;
        winBuffPos = 0;
        winBuffMinPos = -1;
        kmerSpan = 0;
        hpEvents.clear();
        for (int32_t i = 0; i < winSize; ++i) {
            winBuff[i] = Seed();
            winPosSet[i] = false;
        }
        for (size_t i = 0; i < seqBuffer.size(); ++i) {
            seqBuffer[i] = 0x0;
            seqBufferRC[i] = 0x0;
            numBasesIn[i] = 0;
        }
    }
};

int GenerateMinimizers(std::vector<PacBio::Pancake::Int128t>& minimizers, const uint8_t* seq,
                       const int32_t seqLen, const int32_t seqOffset, const int32_t seqId,
                       const int32_t kmerSize, const int32_t winSize, const int32_t spacing,
                       const bool useReverseComplement, const bool useHPC)
{
    const int32_t MAX_SEED_SPAN = 255;

    // Sanity check that the seq is not NULL;
    if (!seq) {
        throw std::runtime_error("Cannot generate minimizers. The sequence is NULL.");
    }
    // Sanity check for the size of kmer. It can't
    // be too large, or it won't fit the uint64_t buffer.
    if (kmerSize <= 0 || kmerSize > 28) {
        throw std::runtime_error(
            "Cannot generate minimizers. The kmerSize is out of bounds, should be in range [0, "
            "28]. kmerSize = " +
            std::to_string(kmerSize));
    }
    if (winSize > MAX_WINDOW_BUFFER_SIZE) {
        throw std::runtime_error(
            "Cannot generate minimizers. The winSize is out of bounds, should be in range [0, " +
            std::to_string(MAX_WINDOW_BUFFER_SIZE) + "]. winSize = " + std::to_string(winSize));
    }
    if (spacing > MAX_SPACING_IN_SEED) {
        throw std::runtime_error(
            "Cannot generate minimizers. The spacing is out of bounds, should be in range [0, " +
            std::to_string(MAX_SPACING_IN_SEED) + "]. spacing = " + std::to_string(spacing));
    }
    // Not technically an error if the seqLen is 0,
    // but there's nothing to do, so return.
    if (seqLen == 0) {
        return 0;
    }

    const minkey_t mask =
        (kmerSize < 32) ? ((((uint64_t)1) << (2 * kmerSize)) - 1)
                        : 0xFFFFFFFFFFFFFFFF;  // Mask the number of required bits for the kmer.

    SpacedBuffer bufferWD;
    int32_t minInBases = kmerSize * (spacing + 1) - spacing;
    int32_t seedSpanSize = minInBases;

    for (int32_t pos = 0, space = 0; pos < seqLen; ++pos, ++space) {
        int8_t b = seq[pos];
        Seed newSeed;  // (UINT64_T_MAX, INT32_T_MAX, INT32_T_MAX, INT8_T_MAX);
        bool newSeedSet = false;

        if (space > spacing) {
            space = 0;
        }

        if (IsNucleotide[b]) {
            // If enabled, skip same bases.
            if (useHPC) {
                // In this case, find the stretch of homopolymers (at least 1 base), and
                // add that stretch to "seed_span" (the seed_span is iteratively increased every
                // loop iteration).
                int32_t hpLen = 1;
                for (hpLen = 1; (pos + hpLen) < seqLen; ++hpLen) {
                    if (seq[pos + hpLen] != b) {
                        break;
                    }
                }
                bufferWD.hpEvents.push_back(hpLen);
                pos += hpLen - 1;
                bufferWD.kmerSpan += hpLen;
                // We need to keep track of the length of each HP event. They contribute as a single
                // base, but in reality, the span is larger. This deque keeps track of all hp events
                // in the current kmer.
                while (bufferWD.hpEvents.size() > static_cast<size_t>(kmerSize)) {
                    bufferWD.kmerSpan -= bufferWD.hpEvents.front();
                    bufferWD.hpEvents.pop_front();
                }
            } else {
                bufferWD.kmerSpan = std::min(
                    (bufferWD.numBasesIn[space] + 1) * (1 + spacing) - spacing, seedSpanSize);
            }

            // Add the base to the buffer.
            bufferWD.seqBuffer[space] =
                ((bufferWD.seqBuffer[space] << 2) | ((((uint64_t)BaseToTwobit[b]))));
            bufferWD.seqBufferRC[space] =
                (bufferWD.seqBufferRC[space] >> 2) |
                ((((uint64_t)BaseToTwobitComplement[b])) << (kmerSize * 2 - 2));
            // Calculate the seed key.
            minkey_t key = bufferWD.seqBuffer[space] & mask;
            minkey_t keyRev = bufferWD.seqBufferRC[space] & mask;

            // Skip symmetric kmers.
            if (bufferWD.seqBuffer[space] == bufferWD.seqBufferRC[space]) {
                continue;
            }

            // Determine the orientation of the key.
            bool isRev = false;
            if (useReverseComplement && keyRev < key) {
                std::swap(key, keyRev);
                isRev = true;
            }

            ++bufferWD.numBasesIn[space];

            if (bufferWD.numBasesIn[space] >= kmerSize && bufferWD.kmerSpan <= MAX_SEED_SPAN) {
                // Minimap2 has another condition here: "kmerSpan < 256". That's because it encodes the kmer span into the seed definition as 8 bits.
                // The 'pos' is the current position which is inclusive. We need to add a +1 to make it non-inclusive, so that the start position is calculated properly.
                int32_t kmerStart = (pos + 1) - (bufferWD.kmerSpan);
                // const int32_t seedSize = bufferWD.kmerSpan * (spacing + 1) - spacing;
                newSeed = Seed(InvertibleHash(key, mask), bufferWD.kmerSpan, seqId,
                               kmerStart + seqOffset, isRev);
                newSeedSet = true;
            }

        } else {
            // If we encountered a non-nucleotide base, this is enough to reset the minimizer coding.
            // First, write out the current minimizer so it doesn't get lost in the void.
            // We only need to write one and not loop through all of them, because when we encounter a seed
            // with the same key then we write out the previous minimizer out. This means that we only need
            // to write out the current minimizer.
            if (bufferWD.numBasesIn[space] >= (winSize + kmerSize) && bufferWD.winBuffMinPos >= 0 &&
                bufferWD.winPosSet[bufferWD.winBuffMinPos]) {
                minimizers.emplace_back(bufferWD.winBuff[bufferWD.winBuffMinPos].To128t());
            }
            bufferWD.Clear(winSize);
        }

        // The first time the buffer is filled, find and add the previous
        // distinct minimizer matches.
        if (bufferWD.numBasesIn[space] == (winSize + kmerSize - 1) && bufferWD.winBuffMinPos >= 0) {
            const auto& minSeed = bufferWD.winBuff[bufferWD.winBuffMinPos];
            // First portion of the circular buffer.
            for (int32_t j = (bufferWD.winBuffPos + 1); j < winSize; ++j) {
                if (bufferWD.winPosSet[j] && minSeed.Compare(bufferWD.winBuff[j].key) == 2) {
                    minimizers.emplace_back(bufferWD.winBuff[j].To128t());
                }
            }
            // Second portion of the circular buffer.
            for (int32_t j = 0; j < bufferWD.winBuffPos; ++j) {
                if (bufferWD.winPosSet[j] && minSeed.Compare(bufferWD.winBuff[j].key) == 2) {
                    minimizers.emplace_back(bufferWD.winBuff[j].To128t());
                }
            }
        }

        if (newSeedSet && bufferWD.winBuffMinPos < 0) {
            // No minimum has been set yet. Set the current buffer pos.
            bufferWD.winBuffMinPos = bufferWD.winBuffPos;
        } else if (newSeedSet && newSeed.key <= bufferWD.winBuff[bufferWD.winBuffMinPos].key) {
            // In this case, even if we encountered the same minimal key, we will write
            // out the previous occurence of this key, and then set the minimum to the current
            // one. This ensures that all equal minimizers in a window are written.
            if (bufferWD.winPosSet[bufferWD.winBuffMinPos] &&
                bufferWD.numBasesIn[space] >= (winSize + kmerSize)) {
                minimizers.emplace_back(bufferWD.winBuff[bufferWD.winBuffMinPos].To128t());
            }
            bufferWD.winBuffMinPos = bufferWD.winBuffPos;
        } else if (bufferWD.winBuffPos == bufferWD.winBuffMinPos &&
                   bufferWD.winBuffMinPos >=
                       0) {  // The entire window has been circled around to the minimum seed key.
            if (bufferWD.winPosSet[bufferWD.winBuffMinPos] &&
                bufferWD.numBasesIn[space] >= (winSize + kmerSize - 1)) {
                minimizers.emplace_back(bufferWD.winBuff[bufferWD.winBuffMinPos].To128t());
            }
            bufferWD.winBuffMinPos = -1;
            // First portion of the circular buffer.
            for (int32_t j = (bufferWD.winBuffPos + 1); j < winSize; ++j) {
                if (bufferWD.winPosSet[j] &&
                    (bufferWD.winBuffMinPos < 0 ||
                     bufferWD.winBuff[bufferWD.winBuffMinPos].key >= bufferWD.winBuff[j].key)) {
                    bufferWD.winBuffMinPos = j;
                }
            }
            // std::cerr << "(1) Looking for new winBuffMinPos, winBuffMinPos = " << winBuffMinPos << "\n";
            // Second portion of the circular buffer.
            for (int32_t j = 0; j < bufferWD.winBuffPos; ++j) {
                if (bufferWD.winPosSet[j] &&
                    (bufferWD.winBuffMinPos < 0 ||
                     bufferWD.winBuff[bufferWD.winBuffMinPos].key >= bufferWD.winBuff[j].key)) {
                    bufferWD.winBuffMinPos = j;
                }
            }
            // std::cerr << "(2) Looking for new winBuffMinPos, winBuffMinPos = " << winBuffMinPos << "\n";
            if (bufferWD.winBuffMinPos < 0 ||
                (newSeedSet && newSeed.key <= bufferWD.winBuff[bufferWD.winBuffMinPos].key)) {
                bufferWD.winBuffMinPos = bufferWD.winBuffPos;
            }
            // std::cerr << "(3) Looking for new winBuffMinPos, winBuffMinPos = " << winBuffMinPos << "\n";
            // Find and add identical seed keys.
            if (bufferWD.numBasesIn[space] >= (winSize + kmerSize - 1)) {
                const auto& minSeed = bufferWD.winBuff[bufferWD.winBuffMinPos];
                // First portion of the circular buffer.
                for (int32_t j = (bufferWD.winBuffPos + 1); j < winSize; ++j) {
                    if (bufferWD.winPosSet[j] && minSeed.Compare(bufferWD.winBuff[j].key) == 2) {
                        minimizers.emplace_back(bufferWD.winBuff[j].To128t());
                    }
                }
                // Second portion of the circular buffer.
                for (int32_t j = 0; j < bufferWD.winBuffPos; ++j) {
                    if (bufferWD.winPosSet[j] && minSeed.Compare(bufferWD.winBuff[j].key) == 2) {
                        minimizers.emplace_back(bufferWD.winBuff[j].To128t());
                    }
                }
            }
        }
        // Store the new seed to the circular buffer.
        bufferWD.winBuff[bufferWD.winBuffPos] = newSeed;
        bufferWD.winPosSet[bufferWD.winBuffPos] = newSeedSet;
        ++bufferWD.winBuffPos;
        if (bufferWD.winBuffPos == winSize) {
            bufferWD.winBuffPos = 0;
        }
    }
    if (bufferWD.winBuffMinPos >= 0 && bufferWD.winPosSet[bufferWD.winBuffMinPos]) {
        minimizers.emplace_back(bufferWD.winBuff[bufferWD.winBuffMinPos].To128t());
    }

    return 0;
}

void GenerateMinimizers(std::vector<PacBio::Pancake::Int128t>& retSeeds,
                        std::vector<int32_t>& retSequenceLengths,
                        const std::vector<FastaSequenceCached>& targetSeqs, const int32_t kmerSize,
                        const int32_t winSize, const int32_t spacing,
                        const bool useReverseComplement, const bool useHPC)
{
    // Collect all seeds for the target sequences.
    retSeeds.clear();
    retSequenceLengths.clear();
    retSequenceLengths.reserve(targetSeqs.size());
    for (int32_t recordId = 0; recordId < static_cast<int32_t>(targetSeqs.size()); ++recordId) {
        const auto& record = targetSeqs[recordId];
        const uint8_t* seq = reinterpret_cast<const uint8_t*>(record.data());
        const int32_t seqId = record.Id();
        int32_t seqLen = record.size();
        retSequenceLengths.emplace_back(seqLen);
        std::vector<PacBio::Pancake::Int128t> newSeeds;
        int rv = GenerateMinimizers(newSeeds, seq, seqLen, 0, seqId, kmerSize, winSize, spacing,
                                    useReverseComplement, useHPC);
        if (rv)
            throw std::runtime_error("Generating minimizers failed for the target sequence, id = " +
                                     std::to_string(recordId));
        retSeeds.insert(retSeeds.end(), newSeeds.begin(), newSeeds.end());
    }
}

void GenerateMinimizers(std::vector<PacBio::Pancake::Int128t>& retSeeds,
                        std::vector<int32_t>& retSequenceLengths,
                        const std::vector<std::string>& targetSeqs, const int32_t kmerSize,
                        const int32_t winSize, const int32_t spacing,
                        const bool useReverseComplement, const bool useHPC)
{
    std::vector<FastaSequenceCached> targetSeqsCached;
    for (int32_t i = 0; i < static_cast<int32_t>(targetSeqs.size()); ++i) {
        targetSeqsCached.emplace_back(
            FastaSequenceCached(std::to_string(i), targetSeqs[i].c_str(), targetSeqs[i].size(), i));
    }
    GenerateMinimizers(retSeeds, retSequenceLengths, targetSeqsCached, kmerSize, winSize, spacing,
                       useReverseComplement, useHPC);
}

}  // namespace SeedDB
}  // namespace Pancake
}  // namespace PacBio
