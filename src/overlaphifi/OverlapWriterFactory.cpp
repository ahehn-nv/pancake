// Authors: Ivan Sovic

#include <pacbio/overlaphifi/OverlapWriterFactory.h>

namespace PacBio {
namespace Pancake {

std::unique_ptr<OverlapWriterBase> OverlapWriterFactory(OverlapWriterType writerType, FILE* fpOut,
                                                        bool writeIds, bool writeCigar)
{
    if (writerType == OverlapWriterType::IPAOvl) {
        return std::unique_ptr<OverlapWriterBase>(
            new OverlapWriterIPAOvl(fpOut, writeIds, writeCigar));
    }
    throw std::runtime_error("Unsupported output format!");
    return nullptr;
}

}  // namespace Pancake
}  // namespace PacBio
