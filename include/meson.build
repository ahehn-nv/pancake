###########
# headers #
###########

if not meson.is_subproject()

  # pacbio/pancake
  install_headers(
    files([
      'pacbio/Version.h',
      'pacbio/pancake/ContiguousFilePart.h',
      'pacbio/alignment/SesDistanceBanded.h',
      'pacbio/overlaphifi/Mapper.h',
      'pacbio/overlaphifi/Overlap.h',
      'pacbio/overlaphifi/OverlapWriter.h',
      'pacbio/overlaphifi/SeedIndex.h',
      'pacbio/seeddb/Minimizers.h',
      'pacbio/seeddb/SeedDBIndexCache.h',
      'pacbio/seeddb/Seed.h',
      'pacbio/seeddb/SeedDBReader.h',
      'pacbio/seeddb/SeedDBReaderCached.h',
      'pacbio/seeddb/SeedDBReaderCachedBlock.h',
      'pacbio/seeddb/SeedDBReaderRawBlock.h',
      'pacbio/seeddb/SeedDBWriter.h',
      'pacbio/seeddb/SequenceSeeds.h',
      'pacbio/seeddb/SequenceSeedsCached.h',
      'pacbio/seqdb/Lookups.h',
      'pacbio/seqdb/Range.h',
      'pacbio/seqdb/FastaSequenceCached.h',
      'pacbio/seqdb/FastaSequenceId.h',
      'pacbio/seqdb/SeqDBIndexCache.h',
      'pacbio/seqdb/SeqDBWriterBase.h',
      'pacbio/seqdb/SeqDBWriter.h',
      'pacbio/seqdb/Twobit.h',
    ])
    # subdir : '')
  )

endif

pancake_include_directories = include_directories(['.', 'pacbio', 'lib/sparsehash'])
