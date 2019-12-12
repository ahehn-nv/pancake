CURRENT_BUILD_DIR?=build
ENABLED_TESTS?=true
export ENABLED_TESTS CURRENT_BUILD_DIR

.PHONY: all build conf unit cram modules

build:
	mkdir -p ${CURRENT_BUILD_DIR} && ninja -C "${CURRENT_BUILD_DIR}" -v
	${CURRENT_BUILD_DIR}/tests/pancake_test

conf:
	rm -rf "${CURRENT_BUILD_DIR}"
	bash -vex scripts/ci/configure_fallback.sh

unit:
	ninja -C "${CURRENT_BUILD_DIR}" -v test

cram: modules
	scripts/cram tests/cram/*.t

all: conf build cram

modules:
	git submodule update --init --recursive
