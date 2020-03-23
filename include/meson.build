###########
# headers #
###########

if not meson.is_subproject()

  #'pacbio/Version.h',

  # pacbio/pancake
  install_headers(
    files([
      'pacbio/pancake/ContiguousFilePart.h',
      ]),
      subdir : 'pacbio/pancake')

  # pacbio/alignment
  install_headers(
    files([
      'pacbio/alignment/SesDistanceBanded.h',
      ]),
      subdir : 'pacbio/alignment')

  # pacbio/overlaphifi
  install_headers(
    files([
      'pacbio/overlaphifi/Mapper.h',
      'pacbio/overlaphifi/Overlap.h',
      'pacbio/overlaphifi/OverlapHifiSettings.h',
      'pacbio/overlaphifi/OverlapWriter.h',
      'pacbio/overlaphifi/SeedHit.h',
      'pacbio/overlaphifi/SeedIndex.h',
      ]),
      subdir : 'pacbio/overlaphifi')

  # pacbio/seeddb
  install_headers(
    files([
      'pacbio/seeddb/Minimizers.h',
      'pacbio/seeddb/Seed.h',
      'pacbio/seeddb/SeedDBIndexCache.h',
      'pacbio/seeddb/SeedDBParameters.h',
      'pacbio/seeddb/SeedDBReader.h',
      'pacbio/seeddb/SeedDBReaderCached.h',
      'pacbio/seeddb/SeedDBReaderCachedBlock.h',
      'pacbio/seeddb/SeedDBReaderRawBlock.h',
      'pacbio/seeddb/SeedDBWriter.h',
      'pacbio/seeddb/SequenceSeeds.h',
      'pacbio/seeddb/SequenceSeedsCached.h',
      ]),
      subdir : 'pacbio/seeddb')

  # pacbio/seqdb
  install_headers(
    files([
      'pacbio/seqdb/CompressedSequence.h',
      'pacbio/seqdb/FastaSequenceCached.h',
      'pacbio/seqdb/FastaSequenceId.h',
      'pacbio/seqdb/Lookups.h',
      'pacbio/seqdb/Range.h',
      'pacbio/seqdb/SeqDBIndexCache.h',
      'pacbio/seqdb/SeqDBReader.h',
      'pacbio/seqdb/SeqDBReaderBase.h',
      'pacbio/seqdb/SeqDBReaderCached.h',
      'pacbio/seqdb/SeqDBReaderCachedBlock.h',
      'pacbio/seqdb/SeqDBWriter.h',
      'pacbio/seqdb/SeqDBWriterBase.h',
      'pacbio/seqdb/Twobit.h',
      'pacbio/seqdb/Util.h',
      ]),
      subdir : 'pacbio/seqdb')

  # pacbio/alignment
  install_headers(
    files([
      'pacbio/util/CommonTypes.h',
      'pacbio/util/TicToc.h',
      ]),
      subdir : 'pacbio/util')

endif

pancake_include_directories = include_directories(['.', 'pacbio', 'lib/sparsehash'])
