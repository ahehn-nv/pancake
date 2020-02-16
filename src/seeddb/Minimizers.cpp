/*
 * Minimizers.cpp
 *
 * This source file was obtained and adjusted from the Raptor
 * graph-based mapping tool codebase.
 *
 *  Created on: Oct 04, 2017
 *      Author: Ivan Sovic
 */

#include <pacbio/seeddb/Minimizers.h>
#include <pacbio/seeddb/Seed.h>
#include <pacbio/seqdb/Lookups.h>
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
            win_pos_set[i] = false;
        }

        for (size_t i = 0; i < seq_buffer.size(); ++i) {
            seq_buffer[i] = 0x0;
            seq_buffer_rc[i] = 0x0;
            num_bases_in[i] = 0;
        }
    }

    ~SpacedBuffer() = default;

    std::array<uint64_t, MAX_SPACING_IN_SEED>
        seq_buffer;  // Holds the current 2-bit seed representation.
    std::array<uint64_t, MAX_SPACING_IN_SEED>
        seq_buffer_rc;  // Holds the reverse complement 2-bit seed at the same position.
    std::array<int32_t, MAX_SPACING_IN_SEED> num_bases_in;  // Number of bases added to the buffer.

    Seed win_buff[MAX_WINDOW_BUFFER_SIZE];  // Define the new circular buffer for the window.
    bool win_pos_set[MAX_WINDOW_BUFFER_SIZE];
    int32_t win_buff_pos = 0;
    int32_t win_buff_min_pos = -1;
    int32_t kmer_span = 0;
    std::deque<int32_t> hp_events;

    void Clear(int32_t winSize)
    {
        // buffer = 0;
        // buffer_rc = 0;
        // num_bases_in = 0;
        win_buff_pos = 0;
        win_buff_min_pos = -1;
        kmer_span = 0;
        hp_events.clear();
        for (int32_t i = 0; i < winSize; ++i) {
            win_buff[i] = Seed();
            win_pos_set[i] = false;
        }
        for (size_t i = 0; i < seq_buffer.size(); ++i) {
            seq_buffer[i] = 0x0;
            seq_buffer_rc[i] = 0x0;
            num_bases_in[i] = 0;
        }
    }
};

int GenerateMinimizers(std::vector<__int128>& minimizers, const uint8_t* seq, int32_t seqLen,
                       int32_t seqOffset, int32_t seqId, int32_t kmerSize, int32_t winSize,
                       int32_t spacing, bool useRC, bool useHPC, int32_t maxHPCLen)
{

    // Sanity check that the seq is not NULL;
    if (!seq) {
        return 1;
    }
    // Sanity check for the size of kmer. It can't
    // be too large, or it won't fit the uint64_t buffer.
    if (kmerSize <= 0 || kmerSize > 32) {
        return 2;
    }
    // Not technically an error if the seqLen is 0,
    // but there's nothing to do, so return.
    if (seqLen == 0) {
        return 0;
    }

    const minkey_t mask =
        (kmerSize < 32) ? ((((uint64_t)1) << (2 * kmerSize)) - 1)
                        : 0xFFFFFFFFFFFFFFFF;  // Mask the number of required bits for the kmer.

    SpacedBuffer buffer_wd;
    int32_t minInBases = kmerSize * (spacing + 1) - spacing;
    int32_t seedSpanSize = minInBases;

    for (int32_t pos = 0, space = 0; pos < seqLen; ++pos, ++space) {
        int8_t b = seq[pos];
        Seed new_seed;  // (UINT64_T_MAX, INT32_T_MAX, INT32_T_MAX, INT8_T_MAX);
        bool new_seed_set = false;

        if (space > spacing) {
            space = 0;
        }

        if (IsNucleotide[b]) {
            // If enabled, skip same bases.
            if (useHPC) {
                // In this case, find the stretch of homopolymers (at least 1 base), and
                // add that stretch to "seed_span" (the seed_span is iteratively increased every
                // loop iteration).
                int32_t hp_len = 1;
                for (hp_len = 1; hp_len < maxHPCLen && (pos + hp_len) < seqLen; ++hp_len) {
                    if (seq[pos + hp_len] != b) {
                        break;
                    }
                }
                buffer_wd.hp_events.push_back(hp_len);
                pos += hp_len - 1;
                buffer_wd.kmer_span += hp_len;
                // We need to keep track of the length of each HP event. They contribute as a single
                // base, but in reality, the span is larger. This deque keeps track of all hp events
                // in the current kmer.
                while (buffer_wd.hp_events.size() > static_cast<size_t>(kmerSize)) {
                    buffer_wd.kmer_span -= buffer_wd.hp_events.front();
                    buffer_wd.hp_events.pop_front();
                }
                if (hp_len >= maxHPCLen) {
                    continue;
                }
            } else {
                buffer_wd.kmer_span = std::min(
                    (buffer_wd.num_bases_in[space] + 1) * (1 + spacing) - spacing, seedSpanSize);
            }

            // Add the base to the buffer.
            buffer_wd.seq_buffer[space] =
                ((buffer_wd.seq_buffer[space] << 2) | ((((uint64_t)BaseToTwobit[b]))));
            buffer_wd.seq_buffer_rc[space] =
                (buffer_wd.seq_buffer_rc[space] >> 2) |
                ((((uint64_t)BaseToTwobitComplement[b])) << (kmerSize * 2 - 2));
            // Calculate the seed key.
            minkey_t key = buffer_wd.seq_buffer[space] & mask;
            minkey_t rev_key = buffer_wd.seq_buffer_rc[space] & mask;

            // Skip symmetric kmers.
            if (buffer_wd.seq_buffer[space] == buffer_wd.seq_buffer_rc[space]) {
                continue;
            }

            // Determine the orientation of the key.
            int8_t flag = MINIMIZER_FLAG_DEFAULT_FWD;
            if (useRC && rev_key < key) {
                std::swap(key, rev_key);
                flag = MINIMIZER_FLAG_IS_REV;
            }

            ++buffer_wd.num_bases_in[space];

            if (buffer_wd.num_bases_in[space] >= kmerSize) {
                // Minimap2 has another condition here: "kmer_span < 256". That's because it encodes the kmer span into the seed definition as 8 bits.
                // The 'pos' is the current position which is inclusive. We need to add a +1 to make it non-inclusive, so that the start position is calculated properly.
                int32_t kmer_start = (pos + 1) - (buffer_wd.kmer_span);
                new_seed = Seed(key, seqId, kmer_start + seqOffset, flag);
                new_seed_set = true;
            }

        } else {
            // If we encountered a non-nucleotide base, this is enough to reset the minimizer coding.
            // First, write out the current minimizer so it doesn't get lost in the void.
            // We only need to write one and not loop through all of them, because when we encounter a seed
            // with the same key then we write out the previous minimizer out. This means that we only need
            // to write out the current minimizer.
            if (buffer_wd.num_bases_in[space] >= (winSize + kmerSize) &&
                buffer_wd.win_buff_min_pos >= 0) {
                minimizers.emplace_back(buffer_wd.win_buff[buffer_wd.win_buff_min_pos].To128t());
            }
            buffer_wd.Clear(winSize);
        }

        // The first time the buffer is filled, find and add the previous
        // distinct minimizer matches.
        if (buffer_wd.num_bases_in[space] == (winSize + kmerSize - 1) &&
            buffer_wd.win_buff_min_pos >= 0) {
            const auto& min_seed = buffer_wd.win_buff[buffer_wd.win_buff_min_pos];
            // First portion of the circular buffer.
            for (int32_t j = (buffer_wd.win_buff_pos + 1); j < winSize; ++j) {
                if (buffer_wd.win_pos_set[j] && min_seed.Compare(buffer_wd.win_buff[j].key) == 2) {
                    minimizers.emplace_back(buffer_wd.win_buff[j].To128t());
                }
            }
            // Second portion of the circular buffer.
            for (int32_t j = 0; j < buffer_wd.win_buff_pos; ++j) {
                if (buffer_wd.win_pos_set[j] && min_seed.Compare(buffer_wd.win_buff[j].key) == 2) {
                    minimizers.emplace_back(buffer_wd.win_buff[j].To128t());
                }
            }
        }

        if (new_seed_set && buffer_wd.win_buff_min_pos < 0) {
            // No minimum has been set yet. Set the current buffer pos.
            buffer_wd.win_buff_min_pos = buffer_wd.win_buff_pos;
        } else if (new_seed_set &&
                   new_seed.key <= buffer_wd.win_buff[buffer_wd.win_buff_min_pos].key) {
            // In this case, even if we encountered the same minimal key, we will write
            // out the previous occurence of this key, and then set the minimum to the current
            // one. This ensures that all equal minimizers in a window are written.
            if (buffer_wd.num_bases_in[space] >= (winSize + kmerSize)) {
                minimizers.emplace_back(buffer_wd.win_buff[buffer_wd.win_buff_min_pos].To128t());
            }
            buffer_wd.win_buff_min_pos = buffer_wd.win_buff_pos;
        } else if (buffer_wd.win_buff_pos == buffer_wd.win_buff_min_pos &&
                   buffer_wd.win_buff_min_pos >=
                       0) {  // The entire window has been circled around to the minimum seed key.
            if (buffer_wd.num_bases_in[space] >= (winSize + kmerSize - 1)) {
                minimizers.emplace_back(buffer_wd.win_buff[buffer_wd.win_buff_min_pos].To128t());
            }
            buffer_wd.win_buff_min_pos = -1;
            // First portion of the circular buffer.
            for (int32_t j = (buffer_wd.win_buff_pos + 1); j < winSize; ++j) {
                if (buffer_wd.win_pos_set[j] &&
                    (buffer_wd.win_buff_min_pos < 0 ||
                     buffer_wd.win_buff[buffer_wd.win_buff_min_pos].key >=
                         buffer_wd.win_buff[j].key)) {
                    buffer_wd.win_buff_min_pos = j;
                }
            }
            // std::cerr << "(1) Looking for new win_buff_min_pos, win_buff_min_pos = " << win_buff_min_pos << "\n";
            // Second portion of the circular buffer.
            for (int32_t j = 0; j < buffer_wd.win_buff_pos; ++j) {
                if (buffer_wd.win_pos_set[j] &&
                    (buffer_wd.win_buff_min_pos < 0 ||
                     buffer_wd.win_buff[buffer_wd.win_buff_min_pos].key >=
                         buffer_wd.win_buff[j].key)) {
                    buffer_wd.win_buff_min_pos = j;
                }
            }
            // std::cerr << "(2) Looking for new win_buff_min_pos, win_buff_min_pos = " << win_buff_min_pos << "\n";
            if (buffer_wd.win_buff_min_pos < 0 ||
                (new_seed_set &&
                 new_seed.key <= buffer_wd.win_buff[buffer_wd.win_buff_min_pos].key)) {
                buffer_wd.win_buff_min_pos = buffer_wd.win_buff_pos;
            }
            // std::cerr << "(3) Looking for new win_buff_min_pos, win_buff_min_pos = " << win_buff_min_pos << "\n";
            // Find and add identical seed keys.
            if (buffer_wd.num_bases_in[space] >= (winSize + kmerSize - 1)) {
                const auto& min_seed = buffer_wd.win_buff[buffer_wd.win_buff_min_pos];
                // First portion of the circular buffer.
                for (int32_t j = (buffer_wd.win_buff_pos + 1); j < winSize; ++j) {
                    if (buffer_wd.win_pos_set[j] &&
                        min_seed.Compare(buffer_wd.win_buff[j].key) == 2) {
                        minimizers.emplace_back(buffer_wd.win_buff[j].To128t());
                    }
                }
                // Second portion of the circular buffer.
                for (int32_t j = 0; j < buffer_wd.win_buff_pos; ++j) {
                    if (buffer_wd.win_pos_set[j] &&
                        min_seed.Compare(buffer_wd.win_buff[j].key) == 2) {
                        minimizers.emplace_back(buffer_wd.win_buff[j].To128t());
                    }
                }
            }
        }
        // Store the new seed to the circular buffer.
        buffer_wd.win_buff[buffer_wd.win_buff_pos] = new_seed;
        buffer_wd.win_pos_set[buffer_wd.win_buff_pos] = new_seed_set;
        ++buffer_wd.win_buff_pos;
        if (buffer_wd.win_buff_pos == winSize) {
            buffer_wd.win_buff_pos = 0;
        }
    }
    if (buffer_wd.win_buff_min_pos >= 0) {
        minimizers.emplace_back(buffer_wd.win_buff[buffer_wd.win_buff_min_pos].To128t());
    }

    return 0;
}

}  // namespace SeedDB
}  // namespace Pancake
}  // namespace PacBio
