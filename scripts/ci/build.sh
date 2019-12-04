#!/usr/bin/env bash
set -vex

#########
# BUILD #
#########

# build
ninja -C "${CURRENT_BUILD_DIR:-build}" -v
