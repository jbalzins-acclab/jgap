#!/usr/bin/env zsh

BUILD_DIR="${BUILD_DIR:-build-debug}"

NINJA_OPT=()
if command -v ninja >/dev/null 2>&1; then
    NINJA_OPT=(-GNinja)
fi

COMPILER_OPTS=()
if [[ -n "${CC:-}" ]]; then
    COMPILER_OPTS+=("-DCMAKE_C_COMPILER=${CC}")
fi
if [[ -n "${CXX:-}" ]]; then
    COMPILER_OPTS+=("-DCMAKE_CXX_COMPILER=${CXX}")
fi

echo "==> Configuring Debug build in '${BUILD_DIR}'..."
cmake -B "${BUILD_DIR}" "${NINJA_OPT[@]}" "${COMPILER_OPTS[@]}" -DCMAKE_BUILD_TYPE=Debug "$@"

echo "==> Building Debug target..."
cmake --build "${BUILD_DIR}" --parallel
echo "==> Successfully built jgap (Debug) in '${BUILD_DIR}'"
