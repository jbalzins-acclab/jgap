#!/bin/bash
set -e

INSTALL_DIR="$(pwd)/quip_install"

# Check if QUIP is already installed
# Find if any quip executable exists in the build directories
if [ -d "$INSTALL_DIR/QUIP/build" ] && find "$INSTALL_DIR/QUIP/build" -name "quip" -type f | read -r; then
    echo "QUIP is already installed at $INSTALL_DIR."
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

echo "Cloning QUIP..."
git clone --recursive https://github.com/libAtoms/QUIP.git

cd QUIP
export QUIP_ROOT="$(pwd)"

echo "Checking out QUIP at a fixed point in time (to freeze the version)..."
git fetch origin
git checkout e5de082e4ba920516271de495b7409208bb8121b
git submodule update --init --recursive

echo "Replacing GAP submodule with the requested branch..."
rm -rf src/GAP
git clone -b tabgap_stuff https://github.com/jesperbygg/GAP.git src/GAP
cd src/GAP
echo "Checking out GAP at a fixed point in time..."
git checkout 33b794686275ff219f7feb2ede962dbff223faff
echo "Initializing GAP submodules..."
git submodule update --init --recursive
cd ../..

echo "Determining architecture..."
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
    export QUIP_ARCH=linux_x86_64_gfortran
    MATH_LIBS="-llapack -lblas"
elif [[ "$OSTYPE" == "darwin"* ]]; then
    export QUIP_ARCH=darwin_x86_64_gfortran
    MATH_LIBS="-framework Accelerate"
else
    echo "Unsupported OS. Defaulting to linux_x86_64_gfortran"
    export QUIP_ARCH=linux_x86_64_gfortran
    MATH_LIBS="-llapack -lblas"
fi

echo "Configuring QUIP for $QUIP_ARCH with fast compilation flags and OpenMP..."
mkdir -p "build/$QUIP_ARCH"

# Create Makefile.inc tailored for fast execution (OpenMP, O3, fast-math)
cat <<EOF > "build/$QUIP_ARCH/Makefile.inc"
# Base architecture settings
include \${QUIP_ROOT}/arch/Makefile.\${QUIP_ARCH}

# Fast compilation and OpenMP
HAVE_GAP=1
HAVE_OPENMP=1
MATH_LINKOPTS=$MATH_LIBS

# Optimization flags
OPT=-O3 -ffast-math -fopenmp
F77_OPT=-O3 -ffast-math -fopenmp -fallow-argument-mismatch
F90_OPT=-O3 -ffast-math -fopenmp -fallow-argument-mismatch
C_OPT=-O3 -ffast-math -fopenmp
CPLUSPLUS_OPT=-O3 -ffast-math -fopenmp
EOF

echo "Running make config to set up SIZEOF_FORTRAN_T and other architecture details..."
# QUIP's make config script requires some user input or default values.
# We'll use yes "" to feed the defaults (empty enters) to the config step.
yes "" | make config

echo "Applying fix for missing compilation rules..."
# Older QUIP Makefile only knows how to compile specifically listed extensions or .f90 files
# in certain directories. GAP's new branch introduces files with .f95 which QUIP's rule
# engine might miss entirely if they are not explicitly handled in QUIP's rules.mk
# We simply rename all .f95 files to .f90 and update the GAP Makefile to expect .o from .f90
find src/GAP -name "*.f95" -exec bash -c 'mv "$0" "${0%.f95}.f90"' {} \;
sed -i.bak 's/\.f95/\.f90/g' src/GAP/Makefile

echo "Building QUIP & GAP..."
make

echo "=================================================="
echo "Installation complete!"
echo "Executables are located in: $QUIP_ROOT/build/$QUIP_ARCH"
echo "You may want to add it to your PATH."
echo "=================================================="
