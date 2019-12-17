###########
# headers #
###########

if not meson.is_subproject()

  # pacbio/pancake
  install_headers(
    files([
      'pacbio/Version.h',
      'pacbio/seqdb/Lookups.h',
      'pacbio/seqdb/Range.h',
      'pacbio/seqdb/SeqDBIndexCache.h',
      'pacbio/seqdb/SeqDBWriterBase.h',
      'pacbio/seqdb/SeqDBWriterCompressed.h',
      'pacbio/seqdb/Twobit.h',
    ])
    # subdir : '')
  )

endif

pancake_include_directories = include_directories(['.', 'pacbio'])