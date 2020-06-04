// Author: Ivan Sovic

#ifndef PANCAKE_OVERLAPHIFI_OVERLAP_WRITER_FACTORY_H
#define PANCAKE_OVERLAPHIFI_OVERLAP_WRITER_FACTORY_H

#include <pacbio/overlaphifi/OverlapWriterBase.h>
#include <pacbio/overlaphifi/OverlapWriterFormat.h>
#include <pacbio/overlaphifi/OverlapWriterIPAOvl.h>
#include <memory>

namespace PacBio {
namespace Pancake {

std::unique_ptr<OverlapWriterBase> OverlapWriterFactory(OverlapWriterFormat writerType, FILE* fpOut,
                                                        bool writeIds, bool writeCigar);
}
}

#endif  // PANCAKE_OVERLAPHIFI_OVERLAP_WRITER_FACTORY_H
