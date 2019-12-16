CURRENT_BUILD_DIR?=build
ENABLED_TESTS?=true
export ENABLED_TESTS CURRENT_BUILD_DIR

.PHONY: all build conf unit cram modules check-formatting

build:
	mkdir -p ${CURRENT_BUILD_DIR} && ninja -C "${CURRENT_BUILD_DIR}" -v test

conf:
	rm -rf "${CURRENT_BUILD_DIR}"
	bash -vex scripts/ci/configure_fallback.sh

unit:
	ninja -C "${CURRENT_BUILD_DIR}" -v test

cram: modules
	scripts/cram tests/cram/*.t

check-formatting:
	tools/check-formatting --all

all: conf build cram check-formatting

modules:
	git submodule update --init --recursive
