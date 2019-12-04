###########
# headers #
###########

if not meson.is_subproject()

  # pacbio/pancake
  install_headers(
    files([
      'pacbio/pancake/Pancake.h',
      'pacbio/pancake/Settings.h',
      'pacbio/pancake/Version.h'
    ]),
    subdir : 'pacbio/pancake')

endif

pancake_include_directories = include_directories('.')