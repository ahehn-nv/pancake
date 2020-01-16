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
#include <deque>

namespace PacBio {
namespace Pancake {
namespace SeedDB {

using minkey_t = uint64_t;

int GenerateMinimizers(std::vector<__int128>& minimizers, const uint8_t* seq, int32_t seqLen,
                       int32_t seqOffset, int32_t seqId, int32_t k, int32_t w, bool useRC,
                       bool useHPC, int32_t maxHPCLen)
{

    // Sanity check that the seq is not NULL;
    if (!seq) {
        return 1;
    }
    // Sanity check for the size of k. It can't
    // be too large, or it won't fit the uint64_t buffer.
    if (k <= 0 || k >= 31) {
        return 2;
    }
    // Not technically an error if the seqLen is 0,
    // but there's nothing to do, so return.
    if (seqLen == 0) {
        return 0;
    }

    const minkey_t mask =
        (((uint64_t)1) << (2 * k)) - 1;  // Mask the number of required bits for the k-mer.
    minkey_t buffer = 0x0;               // Holds the current 2-bit seed representation.
    minkey_t buffer_rc = 0x0;  // Holds the reverse complement 2-bit seed at the same position.
    int32_t num_bases_in = 0;  // Number of bases added to the buffer.
    Seed win_buff[512];        // Define the new circular buffer for the window.
    bool win_pos_set[512];
    int32_t win_buff_pos = 0;
    int32_t win_buff_min_pos = -1;
    int32_t kmer_span = 0;
    std::deque<int32_t> hp_events;

    for (size_t i = 0; i < 512; ++i) {
        win_pos_set[i] = false;
    }

    for (int32_t pos = 0; pos < seqLen; ++pos) {
        int8_t b = seq[pos];
        Seed new_seed;  // (UINT64_T_MAX, INT32_T_MAX, INT32_T_MAX, INT8_T_MAX);
        bool new_seed_set = false;

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
                hp_events.push_back(hp_len);
                pos += hp_len - 1;
                kmer_span += hp_len;
                // We need to keep track of the length of each HP event. They contribute as a single
                // base, but in reality, the span is larger. This deque keeps track of all hp events
                // in the current k-mer.
                while (hp_events.size() > static_cast<size_t>(k)) {
                    kmer_span -= hp_events.front();
                    hp_events.pop_front();
                }
                if (hp_len >= maxHPCLen) {
                    continue;
                }
            } else {
                kmer_span = std::min(num_bases_in + 1, k);
            }

            // Add the base to the buffer.
            buffer = ((buffer << 2) | ((((uint64_t)BaseToTwobit[b]))));
            buffer_rc = (buffer_rc >> 2) | ((((uint64_t)BaseToTwobitComplement[b])) << (k * 2 - 2));
            // Calculate the seed key.
            minkey_t key = buffer & mask;
            minkey_t rev_key = buffer_rc & mask;
            // Skip symmetric k-mers.
            if (buffer == buffer_rc) {
                continue;
            }

            // Determine the orientation of the key.
            int8_t flag = MINIMIZER_FLAG_DEFAULT_FWD;
            if (useRC && rev_key < key) {
                std::swap(key, rev_key);
                flag = MINIMIZER_FLAG_IS_REV;
            }

            ++num_bases_in;

            if (num_bases_in >=
                k) {  // Minimap2 has another condition here: "kmer_span < 256". That's because it encodes the kmer span into the seed definition as 8 bits.
                int32_t kmer_start =
                    (pos + 1) -
                    kmer_span;  // The 'pos' is the current position which is inclusive. We need to add a +1 to make it non-inclusive, so that the start position is calculated properly.
                new_seed = Seed(key, seqId, kmer_start + seqOffset, flag);
                new_seed_set = true;
            }
        } else {
            // If we encountered a non-nucleotide base, this is enough to reset the minimizer coding.
            // First, write out the current minimizer so it doesn't get lost in the void.
            // We only need to write one and not loop through all of them, because when we encounter a seed
            // with the same key then we write out the previous minimizer out. This means that we only need
            // to write out the current minimizer.
            if (num_bases_in >= (w + k) && win_buff_min_pos >= 0) {
                minimizers.emplace_back(win_buff[win_buff_min_pos].To128t());
            }
            num_bases_in = 0;
            win_buff_pos = 0;
            win_buff_min_pos = -1;
            buffer = 0;
            buffer_rc = 0;
            hp_events.clear();
            for (int32_t i = 0; i < w; ++i) {
                win_pos_set[i] = false;
                win_buff[i] = Seed();
            }
        }

        // The first time the buffer is filled, find and add the previous
        // distinct minimizer matches.
        if (num_bases_in == (w + k - 1) && win_buff_min_pos >= 0) {
            const auto& min_seed = win_buff[win_buff_min_pos];
            // First portion of the circular buffer.
            for (int32_t j = (win_buff_pos + 1); j < w; ++j) {
                if (win_pos_set[j] && min_seed.Compare(win_buff[j].key) == 2) {
                    minimizers.emplace_back(win_buff[j].To128t());
                }
            }
            // Second portion of the circular buffer.
            for (int32_t j = 0; j < win_buff_pos; ++j) {
                if (win_pos_set[j] && min_seed.Compare(win_buff[j].key) == 2) {
                    minimizers.emplace_back(win_buff[j].To128t());
                }
            }
        }

        if (new_seed_set && win_buff_min_pos < 0) {
            // No minimum has been set yet. Set the current buffer pos.
            win_buff_min_pos = win_buff_pos;
        } else if (new_seed_set && new_seed.key <= win_buff[win_buff_min_pos].key) {
            // In this case, even if we encountered the same minimal key, we will write
            // out the previous occurence of this key, and then set the minimum to the current
            // one. This ensures that all equal minimizers in a window are written.
            if (num_bases_in >= (w + k)) {
                minimizers.emplace_back(win_buff[win_buff_min_pos].To128t());
            }
            win_buff_min_pos = win_buff_pos;
        } else if (win_buff_pos == win_buff_min_pos &&
                   win_buff_min_pos >=
                       0) {  // The entire window has been circled around to the minimum seed key.
            if (num_bases_in >= (w + k - 1)) {
                minimizers.emplace_back(win_buff[win_buff_min_pos].To128t());
            }
            win_buff_min_pos = -1;
            // First portion of the circular buffer.
            for (int32_t j = (win_buff_pos + 1); j < w; ++j) {
                if (win_pos_set[j] &&
                    (win_buff_min_pos < 0 || win_buff[win_buff_min_pos].key >= win_buff[j].key)) {
                    win_buff_min_pos = j;
                }
            }
            // std::cerr << "(1) Looking for new win_buff_min_pos, win_buff_min_pos = " << win_buff_min_pos << "\n";
            // Second portion of the circular buffer.
            for (int32_t j = 0; j < win_buff_pos; ++j) {
                if (win_pos_set[j] &&
                    (win_buff_min_pos < 0 || win_buff[win_buff_min_pos].key >= win_buff[j].key)) {
                    win_buff_min_pos = j;
                }
            }
            // std::cerr << "(2) Looking for new win_buff_min_pos, win_buff_min_pos = " << win_buff_min_pos << "\n";
            if (win_buff_min_pos < 0 ||
                (new_seed_set && new_seed.key <= win_buff[win_buff_min_pos].key)) {
                win_buff_min_pos = win_buff_pos;
            }
            // std::cerr << "(3) Looking for new win_buff_min_pos, win_buff_min_pos = " << win_buff_min_pos << "\n";
            // Find and add identical seed keys.
            if (num_bases_in >= (w + k - 1)) {
                const auto& min_seed = win_buff[win_buff_min_pos];
                // First portion of the circular buffer.
                for (int32_t j = (win_buff_pos + 1); j < w; ++j) {
                    if (win_pos_set[j] && min_seed.Compare(win_buff[j].key) == 2) {
                        minimizers.emplace_back(win_buff[j].To128t());
                    }
                }
                // Second portion of the circular buffer.
                for (int32_t j = 0; j < win_buff_pos; ++j) {
                    if (win_pos_set[j] && min_seed.Compare(win_buff[j].key) == 2) {
                        minimizers.emplace_back(win_buff[j].To128t());
                    }
                }
            }
        }
        // Store the new seed to the circular buffer.
        win_buff[win_buff_pos] = new_seed;
        win_pos_set[win_buff_pos] = new_seed_set;
        ++win_buff_pos;
        if (win_buff_pos == w) {
            win_buff_pos = 0;
        }
    }
    if (win_buff_min_pos >= 0) {
        minimizers.emplace_back(win_buff[win_buff_min_pos].To128t());
    }

    return 0;
}

}  // namespace SeedDB
}  // namespace Pancake
}  // namespace PacBio
