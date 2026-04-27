#!/bin/bash
set -e

echo "Starting jgap compilation setup..."

echo "=================================================="
echo "Dependency Management"
echo "=================================================="
echo "If you already have all dependencies installed (Eigen3, nlohmann_json, TBB, HighFive, pugixml, OpenBLAS),"
echo "you can simply run: cmake -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j"
echo "If you wish to manage dependencies manually, you can skip vcpkg by running this script with --no-vcpkg"
echo ""
echo "Example of installing dependencies via apt (Ubuntu/Debian):"
echo "   sudo apt install nlohmann-json3-dev libtbb-dev libhdf5-dev libpugixml-dev libopenblas-dev"
echo "   (Note: Eigen3 and HighFive might need to be downloaded manually depending on repository versions)"
echo ""
echo "Example of installing dependencies via Homebrew (macOS):"
echo "   brew install eigen nlohmann-json tbb hdf5 pugixml openblas"
echo "   (Note: HighFive might need to be downloaded manually)"
echo "=================================================="

# Check if vcpkg exists either locally or through environment variable
VCPKG_EXEC=""
USE_VCPKG="true"

# Parse arguments to allow skipping vcpkg
if [ "$1" = "--no-vcpkg" ]; then
    echo "Skipping vcpkg as requested by --no-vcpkg flag."
    USE_VCPKG="false"
fi

if [ "$USE_VCPKG" = "true" ]; then
    if [ -n "$VCPKG_ROOT" ] && [ -f "$VCPKG_ROOT/vcpkg" ]; then
        echo "Found vcpkg via VCPKG_ROOT: $VCPKG_ROOT"
        VCPKG_EXEC="$VCPKG_ROOT/vcpkg"
        LOCAL_VCPKG_ROOT="$VCPKG_ROOT"
    elif command -v vcpkg >/dev/null 2>&1; then
        echo "Found vcpkg in PATH"
        VCPKG_EXEC="vcpkg"
        LOCAL_VCPKG_ROOT=$(dirname $(command -v vcpkg))
    elif [ -d "deps_installed" ] && [ -f "deps_installed/vcpkg" ]; then
        echo "Found local deps_installed/vcpkg folder."
        VCPKG_EXEC="./deps_installed/vcpkg"
        LOCAL_VCPKG_ROOT="$(pwd)/deps_installed"
    else
        echo "vcpkg not found. Installing vcpkg locally into 'deps_installed'..."
        git clone https://github.com/microsoft/vcpkg.git deps_installed
        cd deps_installed
        ./bootstrap-vcpkg.sh -disableMetrics
        cd ..
        VCPKG_EXEC="./deps_installed/vcpkg"
        LOCAL_VCPKG_ROOT="$(pwd)/deps_installed"
    fi

    export VCPKG_ROOT="$LOCAL_VCPKG_ROOT"
    echo "Using VCPKG_ROOT = $VCPKG_ROOT"
    CMAKE_EXTRA_ARGS="-DCMAKE_TOOLCHAIN_FILE=$VCPKG_ROOT/scripts/buildsystems/vcpkg.cmake"
else
    CMAKE_EXTRA_ARGS=""
fi

echo "Configuring CMake for jgap Release build..."
mkdir -p build
cd build

cmake .. \
    $CMAKE_EXTRA_ARGS \
    -DCMAKE_BUILD_TYPE=Release

echo "Compiling jgap..."
# Determine number of CPU cores for parallel build
if command -v nproc > /dev/null; then
    CORES=$(nproc)
elif command -v sysctl > /dev/null; then
    CORES=$(sysctl -n hw.ncpu)
else
    CORES=4
fi

make -j"$CORES"

echo "=================================================="
echo "jgap compiled successfully in Release mode!"
echo "Executables are located in the build/ directory."
echo "=================================================="
