// Authors: Ivan Sovic

#include <pacbio/util/RunLengthEncoding.h>

namespace PacBio {
namespace Pancake {

void RunLengthEncoding(const std::string& seq, std::string& encodedSeq,
                       std::vector<int32_t>& runLengths)
{
    RunLengthEncoding(seq.c_str(), seq.size(), encodedSeq, runLengths);
}

void RunLengthEncoding(const char* seq, int64_t seqLen, std::string& encodedSeq,
                       std::vector<int32_t>& runLengths)
{
    encodedSeq = std::string(seq, seqLen);
    runLengths.clear();
    if (seq == NULL || seqLen == 0) return;

    int64_t outSeqLen = RunLengthEncoding(&encodedSeq[0], seqLen, runLengths);
    encodedSeq[outSeqLen] = '\n';
    encodedSeq.resize(outSeqLen);
}

int64_t RunLengthEncoding(char* seq, int64_t seqLen, std::vector<int32_t>& runLengths)
{
    runLengths.clear();
    if (seq == NULL || seqLen == 0) return 0;

    runLengths.resize(seqLen);
    int64_t prev = 0;
    int64_t count = 0;

    for (int64_t i = 0; i < seqLen; ++i) {
        if (seq[prev] != seq[i]) {
            ++prev;
            count = 0;
            seq[prev] = seq[i];
        }
        ++count;
        runLengths[prev] = count;
    }

    ++prev;
    runLengths.resize(prev);

    return prev;
}

}  // namespace Pancake
}  // namespace PacBio
