###########
# headers #
###########

if not meson.is_subproject()

  # pacbio/pancake
  install_headers(
    files([
      'pacbio/Version.h',
      'pacbio/seeddb/Minimizers.h',
      'pacbio/seeddb/SeedDBIndexCache.h',
      'pacbio/seeddb/Seed.h',
      'pacbio/seeddb/SeedDBReaderCached.h',
      'pacbio/seeddb/SeedDBWriter.h',
      'pacbio/seqdb/Lookups.h',
      'pacbio/seqdb/Range.h',
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
