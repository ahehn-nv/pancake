
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

pancake_cpp_sources = files([
    'overlaphifi/OverlapHifiSettings.cpp',
    'overlaphifi/OverlapHifiWorkflow.cpp',
    'seqdb/CompressedSequence.cpp',
    'seqdb/SeqDBSettings.cpp',
    'seqdb/SeqDBWorkflow.cpp',
    'seqdb/SeqDBWriterCompressed.cpp',
    'seqdb/Twobit.cpp',
])

pancake_cpp_sources += pancake_version_sources

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
  install : pancake_lib_install,
  dependencies : pancake_lib_deps,
  include_directories : pancake_include_directories,
  cpp_args : pancake_warning_flags)

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
    cpp_args : pancake_warning_flags)
endif
