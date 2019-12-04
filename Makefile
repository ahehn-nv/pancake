CURRENT_BUILD_DIR?=build
ENABLED_TESTS?=true
export ENABLED_TESTS CURRENT_BUILD_DIR

.PHONY: build conf

build:
	mkdir -p ${CURRENT_BUILD_DIR} && ninja -C "${CURRENT_BUILD_DIR}" -v
	${CURRENT_BUILD_DIR}/tests/pancake_test

conf:
	rm -rf "${CURRENT_BUILD_DIR}"
	bash -vex scripts/ci/configure_fallback.sh

all: conf build
