#!/usr/bin/env bash

type module >& /dev/null || . /mnt/software/Modules/current/init/bash

module purge

module load gcc
module load ccache

module load meson
module load ninja

module load boost
module load cram
module load gtest

case "${bamboo_planRepository_branchName}" in
  master)
    module load pbcopper/master
    ;;
  *)
    module load pbcopper/develop
    ;;
esac
