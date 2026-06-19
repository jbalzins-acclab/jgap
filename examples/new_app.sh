#!/usr/bin/env bash
#
# Scaffolds a new standalone jgap app: creates <AppName>/<AppName>.cpp and a CMakeLists.txt that links
# jgap::jgap, then prints build/run instructions. The CMakeLists template is embedded below (this script
# replaces the old standalone examples/CMakeLists.txt).
#
# It asks for the jgap build directory and whether jgap was built with vcpkg (so it can locate the
# dependency CMake configs). jgap exports its build tree, so no install step is required.

set -euo pipefail

read -rp "App name (e.g. MyApp): " APP
if [[ -z "$APP" ]]; then
    echo "App name must not be empty." >&2
    exit 1
fi

read -rp "jgap build directory: " JGAP_BUILD
if [[ ! -f "$JGAP_BUILD/jgapConfig.cmake" ]]; then
    echo "warning: $JGAP_BUILD/jgapConfig.cmake not found — did you build jgap there?" >&2
fi

read -rp "Was jgap built with vcpkg? [y/N]: " USED_VCPKG
if [[ "$USED_VCPKG" =~ ^[Yy] ]]; then
    DEPS_PREFIX=$(find "$JGAP_BUILD/vcpkg_installed" -maxdepth 1 -mindepth 1 -type d ! -name vcpkg 2>/dev/null | head -1)
    if [[ -z "$DEPS_PREFIX" ]]; then
        echo "Could not find a vcpkg triplet dir under $JGAP_BUILD/vcpkg_installed." >&2
        read -rp "Dependencies prefix (Eigen/TBB/HighFive configs): " DEPS_PREFIX
    fi
else
    read -rp "Dependencies prefix (Eigen/TBB/HighFive configs): " DEPS_PREFIX
fi

PREFIX_PATH="$JGAP_BUILD;$DEPS_PREFIX"

mkdir -p "$APP"

# ---- source stub ----
cat > "$APP/$APP.cpp" <<EOF
#include <iostream>

#include "core/potentials/Potential.hpp"
#include "core/ValuePtr.hpp"
#include "serialization/SerializationRegistry.hpp"
#include "io/log/CurrentLogger.hpp"

using namespace jgap;

int main(int argc, char** argv) {
    CurrentLogger::initDefault({});

    // TODO: your jgap code here, e.g. load a potential:
    //   auto potential = SerializationRegistry<Potential>::deserialize("pot.h5");

    std::cout << "$APP: hello from jgap" << std::endl;
    return 0;
}
EOF

# ---- CMakeLists template ----
cat > "$APP/CMakeLists.txt" <<EOF
cmake_minimum_required(VERSION 3.11)
project($APP CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Linking jgap::jgap is all that is needed: it brings in jgap's headers and its public dependencies.
find_package(jgap CONFIG REQUIRED)

add_executable($APP $APP.cpp)
target_link_libraries($APP PRIVATE jgap::jgap)
EOF

echo
echo "Created $APP/$APP.cpp and $APP/CMakeLists.txt"
echo
echo "Build:"
echo "  cmake -B $APP/build -S $APP -DCMAKE_PREFIX_PATH=\"$PREFIX_PATH\""
echo "  cmake --build $APP/build"
echo
echo "Run:"
echo "  ./$APP/build/$APP"
