CURRENT_BUILD_DIR?=build
CURRENT_DEBUG_BUILD_DIR?=build-debug
ENABLED_TESTS?=true
export ENABLED_TESTS CURRENT_BUILD_DIR

.PHONY: all build conf conf-debug unit cram modules check-formatting build-debug



################
### Release. ###
################
build:
	mkdir -p ${CURRENT_BUILD_DIR} && ninja -C "${CURRENT_BUILD_DIR}" -v test

conf:
	rm -rf "${CURRENT_BUILD_DIR}"
	bash -vex scripts/ci/configure_fallback.sh

all: conf build cram check-formatting

################
### Debug.   ###
################
build-debug:
	mkdir -p ${CURRENT_DEBUG_BUILD_DIR} && ninja -C "${CURRENT_DEBUG_BUILD_DIR}" -v test

conf-debug:
	rm -rf "${CURRENT_DEBUG_BUILD_DIR}"
	bash -vex scripts/ci/configure_debug_fallback.sh ${CURRENT_DEBUG_BUILD_DIR}

debug: conf-debug build-debug

##############
### Tests. ###
##############
unit:
	ninja -C "${CURRENT_BUILD_DIR}" -v test

cram: modules
	scripts/cram tests/cram/*.t

check-formatting:
	tools/check-formatting --all

##############
### Other. ###
##############
modules:
	git submodule update --init --recursive
