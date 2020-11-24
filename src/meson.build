
################
# version info #
################
pancake_ver_arr = meson.project_version().split('.')
pancake_major_version = pancake_ver_arr[0]
pancake_minor_version = pancake_ver_arr[1]
pancake_patch_version = pancake_ver_arr[2]

pancake_config = configuration_data()
pancake_config.set('PANCAKE_VERSION', meson.project_version())
pancake_config.set('PANCAKE_VERSION_MAJOR', pancake_major_version)
pancake_config.set('PANCAKE_VERSION_MINOR', pancake_minor_version)
pancake_config.set('PANCAKE_VERSION_PATCH', pancake_patch_version)

pancake_version_sources = [
  vcs_tag(
    input : 'PancakeGitHash.cpp.in',
    output : 'PancakeGitHash.cpp',
    command : ['git', 'describe', '--always', '--dirty=*'],
    replace_string : '@PANCAKE_GIT_SHA1@'),
  configure_file(
    input : 'PancakeVersion.cpp.in',
    output : 'PancakeVersion.cpp',
    configuration : pancake_config)]

###########
# library #
###########

ksw2_cpp_simd_sources = files([
    'lib/ksw2/ksw2_extd2_sse.cpp',
    'lib/ksw2/ksw2_exts2_sse.cpp',
    'lib/ksw2/ksw2_extz2_sse.cpp',
])
ksw2_cpp_ll_sse_sources = files([
    'lib/ksw2/ksw2_ll_sse.cpp',
])
ksw2_cpp_dispatch_sources = files([
    'lib/ksw2/ksw2_dispatch.cpp',
])

pancake_cpp_sources = files([
    'lib/ksw2/kalloc.cpp',

    'alignment/AlignmentTools.cpp',
    'alignment/SesDistanceBanded.cpp',
    'main/dbfilter/DBFilterSettings.cpp',
    'main/dbfilter/DBFilterWorkflow.cpp',
    'main/overlaphifi/OverlapHifiSettings.cpp',
    'main/overlaphifi/OverlapHifiWorkflow.cpp',
    'main/seeddb/SeedDBSettings.cpp',
    'main/seeddb/SeedDBWorkflow.cpp',
    'main/seqdb/SeqDBSettings.cpp',
    'main/seqdb/SeqDBWorkflow.cpp',
    'main/seqfetch/SeqFetchSettings.cpp',
    'main/seqfetch/SeqFetchWorkflow.cpp',
    'pancake/AlignerBase.cpp',
    'pancake/AlignerKSW2.cpp',
    'pancake/AlignerEdlib.cpp',
    'pancake/AlignerSES1.cpp',
    'pancake/AlignerSES2.cpp',
    'pancake/AlignerFactory.cpp',
    'pancake/AlignmentSeeded.cpp',
    'pancake/CompressedSequence.cpp',
    'pancake/DPChain.cpp',
    'pancake/FastaSequenceId.cpp',
    'pancake/MapperCLR.cpp',
    'pancake/MapperHiFi.cpp',
    'pancake/Minimizers.cpp',
    'pancake/Overlap.cpp',
    'pancake/OverlapWriterBase.cpp',
    'pancake/OverlapWriterFactory.cpp',
    'pancake/OverlapWriterIPAOvl.cpp',
    'pancake/OverlapWriterM4.cpp',
    'pancake/OverlapWriterPAF.cpp',
    'pancake/OverlapWriterSAM.cpp',
    'pancake/Secondary.cpp',
    'pancake/SeedHit.cpp',
    'pancake/SeedHitWriter.cpp',
    'pancake/SeedDBIndexCache.cpp',
    'pancake/SeedDBReader.cpp',
    'pancake/SeedDBReaderCached.cpp',
    'pancake/SeedDBReaderCachedBlock.cpp',
    'pancake/SeedDBReaderRawBlock.cpp',
    'pancake/SeedDBWriter.cpp',
    'pancake/SeedIndex.cpp',
    'pancake/SeqDBIndexCache.cpp',
    'pancake/SeqDBReader.cpp',
    'pancake/SeqDBReaderCached.cpp',
    'pancake/SeqDBReaderCachedBlock.cpp',
    'pancake/SeqDBWriter.cpp',
    'pancake/SequenceSeeds.cpp',
    'pancake/SequenceSeedsCached.cpp',
    'pancake/Twobit.cpp',
    'util/FileIO.cpp',
    'util/RunLengthEncoding.cpp',
    'util/TicToc.cpp',
])

pancake_cpp_sources += pancake_version_sources


##########
# SIMD   #
##########

simd_mod = import('unstable_simd')

ksw2_flags = ['-DHAVE_KALLOC']

cpp = meson.get_compiler('cpp')
foreach cflag: [
  '-Wc++-compat',
  '-Wno-unused-parameter',
  '-Wno-unused-variable',
  '-Wno-sign-compare',
  '-Wno-unused-but-set-variable']
    if cpp.has_argument(cflag)
      ksw2_flags += cflag
    endif
endforeach

ksw2_simd_flag = []
ksw2_simd_libs = []

ksw2_simd_libs += simd_mod.check(
  'll_sse',
  sse2 : ksw2_cpp_ll_sse_sources,
  include_directories : ksw2_include_directories,
  cpp_args : ksw2_flags,
  compiler : cpp)[0]

if opt_sse41
  # SSE 4.1
  ksw2_simd_flag += ['-DKSW_CPU_DISPATCH']

  ksw2_simd_libs += simd_mod.check(
    'ksw2_sse41',
    sse41 : ksw2_cpp_simd_sources,
    include_directories : ksw2_include_directories,
    cpp_args : [ksw2_flags, ksw2_simd_flag],
    compiler : cpp)[0]

  ksw2_simd_libs += simd_mod.check(
    'ksw2_dispatch',
    sse41 : ksw2_cpp_dispatch_sources,
    include_directories : ksw2_include_directories,
    cpp_args : [ksw2_flags, ksw2_simd_flag],
    compiler : cpp)[0]

  ksw2_simd_flag += ['-DKSW_SSE2_ONLY']

endif

# SSE2
ksw2_simd_libs += simd_mod.check(
  'ksw2_sse2',
  sse2 : ksw2_cpp_simd_sources,
  include_directories : ksw2_include_directories,
  cpp_args : [ksw2_flags, ksw2_simd_flag],
  compiler : cpp)[0]

# this is a crutch to unify all compiled objects into
# one library, in order to avoid tons of libraries
# showing up in the pkg-config file
# https://github.com/mesonbuild/meson/issues/2478#issuecomment-351132872
ksw2_simd_all_objects = []
foreach i : ksw2_simd_libs
  ksw2_simd_all_objects += i.extract_all_objects()
endforeach

###########
# Library #
###########
# install library if
# - either running as a proper project
# - or using shared libraries
pancake_lib_install = (not meson.is_subproject()) or (get_option('default_library') == 'shared')

pancake_exe_dir = meson.current_source_dir()

pancake_lib = library(
  'pancake',
  pancake_cpp_sources,
  soversion : meson.project_version(),
  version : meson.project_version(),
  objects : ksw2_simd_all_objects,
  install : pancake_lib_install,
  dependencies : pancake_lib_deps,
  include_directories : pancake_include_directories,
  cpp_args : pancake_flags + ksw2_flags + ksw2_simd_flag)

##############
# executable #
##############

if not meson.is_subproject()
  pancake_main_exe = executable(
    'pancake', files([
      'main/main.cpp',
    ]),
    install : true,
    dependencies : pancake_lib_deps,
    include_directories : pancake_include_directories,
    link_with : pancake_lib,
    cpp_args : pancake_flags + ksw2_flags + ksw2_simd_flag)
endif
