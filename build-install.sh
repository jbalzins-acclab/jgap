#!/usr/bin/env zsh
set -e

PYTHON_BIN="$(which python3 2>/dev/null || which python 2>/dev/null || echo "python3")"

if [[ -n "${VIRTUAL_ENV}" ]]; then
    DEFAULT_PREFIX="${VIRTUAL_ENV}"
else
    DEFAULT_PREFIX="$HOME/.local"
fi

PREFIX="${PREFIX:-${DEFAULT_PREFIX}}"
BUILD_DIR="${BUILD_DIR:-build}"

# Detect if cached build has a mismatched Python or prefix
if [[ -f "${BUILD_DIR}/CMakeCache.txt" ]]; then
    CACHED_PY=$(grep -E '^_?Python3?_EXECUTABLE' "${BUILD_DIR}/CMakeCache.txt" | head -n1 | cut -d'=' -f2 || true)
    CACHED_PREFIX=$(grep -E '^CMAKE_INSTALL_PREFIX:PATH' "${BUILD_DIR}/CMakeCache.txt" | cut -d'=' -f2 || true)
    if [[ -n "${CACHED_PY}" && "${CACHED_PY}" != "${PYTHON_BIN}" ]] || [[ -n "${CACHED_PREFIX}" && "${CACHED_PREFIX}" != "${PREFIX}" ]]; then
        echo "==> Python interpreter or install prefix changed (${CACHED_PY} -> ${PYTHON_BIN}); cleaning '${BUILD_DIR}'..."
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

echo "==> Configuring Release build in '${BUILD_DIR}' with prefix '${PREFIX}'..."
cmake -B "${BUILD_DIR}" "${NINJA_OPT[@]}" "${COMPILER_OPTS[@]}" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="${PREFIX}" \
    -DPython_EXECUTABLE="${PYTHON_BIN}" \
    -DPython3_EXECUTABLE="${PYTHON_BIN}" \
    "$@"

echo "==> Building Release target..."
cmake --build "${BUILD_DIR}" --parallel

echo "==> Installing to '${PREFIX}'..."
cmake --install "${BUILD_DIR}"
echo "==> Successfully installed jgap (Release) to ${PREFIX}"

# Export variables in current shell session & add to shell rc if installing to ~/.local (only when NOT in a venv)
if [[ "${PREFIX}" == "$HOME/.local"* ]]; then
    OS_NAME="$(uname -s)"
    CURRENT_SHELL="$(basename "${SHELL:-/bin/zsh}")"

    RC_FILE=""
    if [[ "${CURRENT_SHELL}" == "zsh" ]]; then
        RC_FILE="$HOME/.zshrc"
    elif [[ "${CURRENT_SHELL}" == "bash" ]]; then
        if [[ "${OS_NAME}" == "Darwin" ]]; then
            RC_FILE="$HOME/.bash_profile"
        else
            RC_FILE="$HOME/.bashrc"
        fi
    else
        RC_FILE="$HOME/.zshrc"
    fi

    DYLD_VAR="LD_LIBRARY_PATH"
    if [[ "${OS_NAME}" == "Darwin" ]]; then
        DYLD_VAR="DYLD_LIBRARY_PATH"
    fi

    PY_SITE_DIR=$("${PYTHON_BIN}" -c "import sys; print(f'lib/python{sys.version_info.major}.{sys.version_info.minor}/site-packages')")

    # Export into active shell session (when sourced)
    export CPLUS_INCLUDE_PATH="$HOME/.local/include:${CPLUS_INCLUDE_PATH:-}"
    export LIBRARY_PATH="$HOME/.local/lib:${LIBRARY_PATH:-}"
    if [[ -z "${VIRTUAL_ENV}" ]]; then
        export PYTHONPATH="$HOME/.local/${PY_SITE_DIR}:${PYTHONPATH:-}"
    fi
    if [[ "${OS_NAME}" == "Darwin" ]]; then
        export DYLD_LIBRARY_PATH="$HOME/.local/lib:${DYLD_LIBRARY_PATH:-}"
    else
        export LD_LIBRARY_PATH="$HOME/.local/lib:${LD_LIBRARY_PATH:-}"
    fi

    echo "==> Exported CPLUS_INCLUDE_PATH, LIBRARY_PATH, and ${DYLD_VAR} in current shell"

    # Persist in shell rc file with dynamic Python version check and virtualenv guard
    touch "${RC_FILE}"
    HEADER_EXPORT='export CPLUS_INCLUDE_PATH="$HOME/.local/include:${CPLUS_INCLUDE_PATH:-}"'
    LIB_EXPORT='export LIBRARY_PATH="$HOME/.local/lib:${LIBRARY_PATH:-}"'
    # Dynamic PYTHONPATH that evaluates for whichever Python is active, but only when not in a venv
    PYTHON_EXPORT='[[ -z "${VIRTUAL_ENV}" ]] && export PYTHONPATH="$HOME/.local/lib/python$(python3 -c "import sys; print(f\"{sys.version_info.major}.{sys.version_info.minor}\")" 2>/dev/null)/site-packages${PYTHONPATH:+:$PYTHONPATH}"'
    DYLD_EXPORT="export ${DYLD_VAR}=\"\$HOME/.local/lib:\${${DYLD_VAR}:-}\""

    # Clean out any old hardcoded PYTHONPATH that broke version switching
    if [[ "${OS_NAME}" == "Darwin" ]]; then
        sed -i '' '/PYTHONPATH.*\.local\/lib\/python3\.[0-9]*/d' "${RC_FILE}" 2>/dev/null || true
    else
        sed -i '/PYTHONPATH.*\.local\/lib\/python3\.[0-9]*/d' "${RC_FILE}" 2>/dev/null || true
    fi

    ADDED=0
    if ! grep -qs 'CPLUS_INCLUDE_PATH.*\.local/include' "${RC_FILE}"; then
        echo "${HEADER_EXPORT}" >> "${RC_FILE}"
        ADDED=1
    fi

    if ! grep -qs 'LIBRARY_PATH.*\.local/lib' "${RC_FILE}"; then
        echo "${LIB_EXPORT}" >> "${RC_FILE}"
        ADDED=1
    fi

    if ! grep -qs 'PYTHONPATH.*\.local/lib' "${RC_FILE}"; then
        echo "${PYTHON_EXPORT}" >> "${RC_FILE}"
        ADDED=1
    fi

    if ! grep -qs "${DYLD_VAR}.*\.local/lib" "${RC_FILE}"; then
        echo "${DYLD_EXPORT}" >> "${RC_FILE}"
        ADDED=1
    fi

    if [[ ${ADDED} -eq 1 ]]; then
        echo "==> Updated $HOME/.local paths in ${RC_FILE}"
    else
        echo "==> ${RC_FILE} is up-to-date"
    fi
fi
