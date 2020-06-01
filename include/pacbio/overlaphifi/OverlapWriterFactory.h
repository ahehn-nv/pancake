// Author: Ivan Sovic

#ifndef PANCAKE_OVERLAPHIFI_OVERLAP_WRITER_FACTORY_H
#define PANCAKE_OVERLAPHIFI_OVERLAP_WRITER_FACTORY_H

#include <pacbio/overlaphifi/OverlapWriterBase.h>
#include <pacbio/overlaphifi/OverlapWriterIPAOvl.h>
#include <memory>

namespace PacBio {
namespace Pancake {

enum class OverlapWriterType
{
    IPAOvl,
    M4,
    M4Extended,
};

std::unique_ptr<OverlapWriterBase> OverlapWriterFactory(OverlapWriterType writerType, FILE* fpOut,
                                                        bool writeIds, bool writeCigar);
}
}

#endif  // PANCAKE_OVERLAPHIFI_OVERLAP_WRITER_FACTORY_H
