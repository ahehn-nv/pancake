unset PKG_CONFIG_LIBDIR
module purge
module load git
module load gcc
module load ccache
module load meson
module load ninja

module load zlib

# Substitute for static versions of our dependencies.
export PKG_CONFIG_LIBDIR=${PKG_CONFIG_LIBDIR//pkgconfig/pkgconfig-static}

# for racon meson
module load cmake gtest

echo "PKG_CONFIG_LIBDIR=${PKG_CONFIG_LIBDIR}"
