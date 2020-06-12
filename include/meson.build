###########
# headers #
###########

if not meson.is_subproject()

  #'pacbio/Version.h',

  # lib/flat_hash_map
  install_headers(
    files([
      'lib/flat_hash_map/bytell_hash_map.hpp',
      'lib/flat_hash_map/flat_hash_map.hpp',
      'lib/flat_hash_map/unordered_map.hpp',
      ]),
      subdir : 'lib/flat_hash_map')

  # lib/kxsort
  install_headers(
    files([
      'lib/kxsort/kxsort.h',
      ]),
      subdir : 'lib/kxsort')

  # lib/sparsehash/internal
  install_headers(
    files([
      'lib/sparsehash/internal/densehashtable.h',
      'lib/sparsehash/internal/hashtable-common.h',
      'lib/sparsehash/internal/libc_allocator_with_realloc.h',
      'lib/sparsehash/internal/sparseconfig.h',
      'lib/sparsehash/internal/sparsehashtable.h',
      ]),
      subdir : 'lib/sparsehash/internal/')

  # lib/sparsehash/sparsehash/internal
  install_headers(
    files([
      'lib/sparsehash/sparsehash/internal/densehashtable.h',
      'lib/sparsehash/sparsehash/internal/hashtable-common.h',
      'lib/sparsehash/sparsehash/internal/libc_allocator_with_realloc.h',
      'lib/sparsehash/sparsehash/internal/sparseconfig.h',
      'lib/sparsehash/sparsehash/internal/sparsehashtable.h',
      ]),
      subdir : 'lib/sparsehash/sparsehash/internal')

  # lib/sparsehash/sparsehash
  install_headers(
    files([
      'lib/sparsehash/sparsehash/template_util.h',
      'lib/sparsehash/sparsehash/type_traits.h',
      ]),
      subdir : 'lib/sparsehash/sparsehash')

  # lib/sparsehash
  install_headers(
    files([
      'lib/sparsehash/template_util.h',
      'lib/sparsehash/type_traits.h',
      ]),
      subdir : 'lib/sparsehash')

  # pacbio/alignment
  install_headers(
    files([
      'pacbio/alignment/AlignmentTools.h',
      'pacbio/alignment/DiffCounts.h',
      'pacbio/alignment/SesAlignBanded.hpp',
      'pacbio/alignment/SesDistanceBanded.h',
      'pacbio/alignment/Ses2DistanceBanded.hpp',
      'pacbio/alignment/SesOptions.h',
      'pacbio/alignment/SesResults.h',
      ]),
      subdir : 'pacbio/alignment')

  # pacbio/overlaphifi
  install_headers(
    files([
      'pacbio/overlaphifi/Mapper.h',
      'pacbio/overlaphifi/Overlap.h',
      'pacbio/overlaphifi/OverlapHifiSettings.h',
      'pacbio/overlaphifi/OverlapWriterBase.h',
      'pacbio/overlaphifi/OverlapWriterFactory.h',
      'pacbio/overlaphifi/OverlapWriterFormat.h',
      'pacbio/overlaphifi/OverlapWriterIPAOvl.h',
      'pacbio/overlaphifi/OverlapWriterM4.h',
      'pacbio/overlaphifi/OverlapWriterSAM.h',
      'pacbio/overlaphifi/SeedHit.h',
      'pacbio/overlaphifi/SeedIndex.h',
      ]),
      subdir : 'pacbio/overlaphifi')

  # pacbio/pancake
  install_headers(
    files([
      'pacbio/pancake/ContiguousFilePart.h',
      ]),
      subdir : 'pacbio/pancake')

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
      'pacbio/util/FileIO.h',
      'pacbio/util/RunLengthEncoding.h',
      'pacbio/util/TicToc.h',
      ]),
      subdir : 'pacbio/util')

endif

pancake_include_directories = include_directories(['.', 'pacbio', 'lib/sparsehash'])
