#!/usr/bin/env zsh

# Allow globs that do not match any files without raising an error
setopt NULL_GLOB NO_NOMATCH 2>/dev/null || true

# Determine script directory whether executed or sourced
if [[ -n "${ZSH_VERSION}" ]]; then
    SCRIPT_DIR="${${(%):-%x}:A:h}"
elif [[ -n "${BASH_SOURCE[0]}" ]]; then
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
else
    SCRIPT_DIR="$(pwd)"
fi

cd "${SCRIPT_DIR}"

echo "==> Cleaning build directories and generated artifacts..."

# Remove CMake build directories (directories only, never .sh files)
rm -rf build/ build-debug/
find . -maxdepth 1 -type d -name "build-*" -exec rm -rf {} + 2>/dev/null || true

# Remove compiled Python extensions in-tree
rm -f python/jgap/*.so python/jgap/*.dylib python/jgap/*.pyd

# Remove Python bytecode and cache
find . -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true
find . -type f -name "*.pyc" -delete 2>/dev/null || true

# Remove generated CMake files in source root if any
rm -rf CMakeFiles CMakeCache.txt cmake_install.cmake generated

echo "==> Clean complete."
