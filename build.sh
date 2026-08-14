#!/usr/bin/env zsh

BUILD_DIR="${BUILD_DIR:-build}"

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

echo "==> Ensuring Python pybind11 dependency is installed..."
python3 -m pip install pybind11 --break-system-packages 2>/dev/null || python3 -m pip install pybind11 || true

echo "==> Configuring Release build in '${BUILD_DIR}'..."
cmake -B "${BUILD_DIR}" "${NINJA_OPT[@]}" "${COMPILER_OPTS[@]}" -DCMAKE_BUILD_TYPE=Release "$@"

echo "==> Building Release target..."
cmake --build "${BUILD_DIR}" --parallel
echo "==> Successfully built jgap (Release) in '${BUILD_DIR}'"
