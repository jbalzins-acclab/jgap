#!/usr/bin/env zsh
set -e

PYTHON_BIN="$(which python3 2>/dev/null || which python 2>/dev/null || echo "python3")"
BUILD_DIR="${BUILD_DIR:-build}"

# Detect if cached build has a mismatched Python
if [[ -f "${BUILD_DIR}/CMakeCache.txt" ]]; then
    CACHED_PY=$(grep -E '^_?Python3?_EXECUTABLE' "${BUILD_DIR}/CMakeCache.txt" | head -n1 | cut -d'=' -f2 || true)
    if [[ -n "${CACHED_PY}" && "${CACHED_PY}" != "${PYTHON_BIN}" ]]; then
        echo "==> Python interpreter changed (${CACHED_PY} -> ${PYTHON_BIN}); cleaning '${BUILD_DIR}'..."
        rm -rf "${BUILD_DIR}"
    fi
fi

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

echo "==> Using Python: ${PYTHON_BIN} ($(${PYTHON_BIN} --version 2>&1))"
echo "==> Ensuring pybind11 is installed for ${PYTHON_BIN}..."
"${PYTHON_BIN}" -m pip install pybind11 --break-system-packages 2>/dev/null || "${PYTHON_BIN}" -m pip install pybind11 || true

echo "==> Configuring Release build in '${BUILD_DIR}'..."
cmake -B "${BUILD_DIR}" "${NINJA_OPT[@]}" "${COMPILER_OPTS[@]}" \
    -DCMAKE_BUILD_TYPE=Release \
    -DPython_EXECUTABLE="${PYTHON_BIN}" \
    -DPython3_EXECUTABLE="${PYTHON_BIN}" \
    "$@"

echo "==> Building Release target..."
cmake --build "${BUILD_DIR}" --parallel
echo "==> Successfully built jgap (Release) in '${BUILD_DIR}'"
