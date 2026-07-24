#!/usr/bin/env zsh

PREFIX="${PREFIX:-$HOME/.local}"
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

echo "==> Configuring Release build in '${BUILD_DIR}'..."
cmake -B "${BUILD_DIR}" "${NINJA_OPT[@]}" "${COMPILER_OPTS[@]}" -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="${PREFIX}" "$@"

echo "==> Building Release target..."
cmake --build "${BUILD_DIR}" --parallel

echo "==> Installing to '${PREFIX}'..."
cmake --install "${BUILD_DIR}"
echo "==> Successfully installed jgap (Release) to ${PREFIX}"

# Export variables in current shell session & add to shell rc
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

    # Export into active shell session (when sourced)
    export CPLUS_INCLUDE_PATH="$HOME/.local/include:${CPLUS_INCLUDE_PATH:-}"
    export LIBRARY_PATH="$HOME/.local/lib:${LIBRARY_PATH:-}"
    if [[ "${OS_NAME}" == "Darwin" ]]; then
        export DYLD_LIBRARY_PATH="$HOME/.local/lib:${DYLD_LIBRARY_PATH:-}"
    else
        export LD_LIBRARY_PATH="$HOME/.local/lib:${LD_LIBRARY_PATH:-}"
    fi

    echo "==> Exported CPLUS_INCLUDE_PATH, LIBRARY_PATH, and ${DYLD_VAR} in current shell"

    # Persist in shell rc file if not present
    touch "${RC_FILE}"
    HEADER_EXPORT='export CPLUS_INCLUDE_PATH="$HOME/.local/include:${CPLUS_INCLUDE_PATH:-}"'
    LIB_EXPORT='export LIBRARY_PATH="$HOME/.local/lib:${LIBRARY_PATH:-}"'
    DYLD_EXPORT="export ${DYLD_VAR}=\"\$HOME/.local/lib:\${${DYLD_VAR}:-}\""

    ADDED=0
    if ! grep -qs 'CPLUS_INCLUDE_PATH.*\.local/include' "${RC_FILE}"; then
        echo "${HEADER_EXPORT}" >> "${RC_FILE}"
        ADDED=1
    fi

    if ! grep -qs 'LIBRARY_PATH.*\.local/lib' "${RC_FILE}"; then
        echo "${LIB_EXPORT}" >> "${RC_FILE}"
        ADDED=1
    fi

    if ! grep -qs "${DYLD_VAR}.*\.local/lib" "${RC_FILE}"; then
        echo "${DYLD_EXPORT}" >> "${RC_FILE}"
        ADDED=1
    fi

    if [[ ${ADDED} -eq 1 ]]; then
        echo "==> Added $HOME/.local paths to ${RC_FILE}"
    else
        echo "==> ${RC_FILE} already contains $HOME/.local paths"
    fi
fi
