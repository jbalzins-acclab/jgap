#!/bin/bash
set -e

INSTALL_DIR="$(pwd)/lammps_install"

# Check if LAMMPS is already installed
if [ -f "$INSTALL_DIR/lammps/build/lmp" ]; then
    echo "LAMMPS is already installed at $INSTALL_DIR/lammps/build/lmp."
    echo "Skipping installation. If you want to reinstall, remove the $INSTALL_DIR directory first."
    exit 0
fi

# Clean up any partial installations
if [ -d "$INSTALL_DIR" ]; then
    echo "Removing incomplete/old installation at $INSTALL_DIR..."
    rm -rf "$INSTALL_DIR"
fi

echo "Setting up installation directory..."
mkdir -p "$INSTALL_DIR"
cd "$INSTALL_DIR"

echo "Cloning LAMMPS..."
git clone https://github.com/lammps/lammps.git

cd lammps
export LAMMPS_ROOT="$(pwd)"

echo "Checking out LAMMPS to specific commit ee93998fe5b05c75e0640aaa4a3b30316bc959bc..."
git fetch origin
git checkout ee93998fe5b05c75e0640aaa4a3b30316bc959bc

echo "Locating TabGAP installation..."
# TabGAP is now installed using the install_tabgap.sh script in the validation/ directory
TABGAP_DIR="../../tabgap_install/tabgap"

if [ ! -d "$TABGAP_DIR" ]; then
    echo "Error: TabGAP not found at $TABGAP_DIR. Please run validation/install_tabgap.sh first."
    exit 1
fi

echo "Copying TabGAP LAMMPS source files into LAMMPS src directory..."
if [ -d "$TABGAP_DIR/lammps" ]; then
    cp "$TABGAP_DIR/lammps/"*.cpp src/
    cp "$TABGAP_DIR/lammps/"*.h src/
else
    echo "Warning: No 'lammps' directory found in TabGAP repo. Copying all .cpp and .h files as a fallback..."
    find "$TABGAP_DIR" -maxdepth 1 -name "*.cpp" -exec cp {} src/ \; || true
    find "$TABGAP_DIR" -maxdepth 1 -name "*.h" -exec cp {} src/ \; || true
fi

echo "Creating build directory..."
mkdir -p build
cd build

echo "Configuring LAMMPS with CMake..."

CMAKE_ARGS=(
    "-DCMAKE_BUILD_TYPE=Release"
    "-DBUILD_SHARED_LIBS=yes"
    "-DPKG_OPENMP=yes"
    "-DBUILD_OMP=yes"
    "-DBUILD_HDF5=yes"
    "-DPKG_ML-GAP=yes"
    "-DDOWNLOAD_EIGEN3=yes"
    "-DCMAKE_CXX_STANDARD=17"
)

# For macOS, Apple Clang does not ship with OpenMP built-in.
# We can use the Homebrew `libomp` package and explicitly point CMake to it.
if [[ "$OSTYPE" == "darwin"* ]]; then
    if command -v brew > /dev/null && brew --prefix libomp > /dev/null 2>&1; then
        LIBOMP_PREFIX=$(brew --prefix libomp)
        echo "Found libomp via Homebrew at $LIBOMP_PREFIX"
        CMAKE_ARGS+=("-DOpenMP_CXX_FLAGS=-Xpreprocessor -fopenmp -I$LIBOMP_PREFIX/include")
        CMAKE_ARGS+=("-DOpenMP_CXX_LIB_NAMES=omp")
        CMAKE_ARGS+=("-DOpenMP_omp_LIBRARY=$LIBOMP_PREFIX/lib/libomp.dylib")
    else
        echo "Warning: libomp not found via Homebrew. If CMake fails to find OpenMP, run: brew install libomp"
    fi
fi

cmake ../cmake "${CMAKE_ARGS[@]}"

echo "Building LAMMPS..."
if command -v nproc > /dev/null; then
    CORES=$(nproc)
elif command -v sysctl > /dev/null; then
    CORES=$(sysctl -n hw.ncpu)
else
    CORES=4
fi

make -j"$CORES"

echo "=================================================="
echo "Installation complete!"
echo "LAMMPS executable is located in: $LAMMPS_ROOT/build/lmp"
echo "You may want to add it to your PATH."
echo "=================================================="
