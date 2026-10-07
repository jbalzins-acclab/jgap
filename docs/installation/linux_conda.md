# Installation Guide: Linux (No `sudo` — Conda / Mamba / Pixi)

This guide covers setting up `jgap` entirely in user space without administrator (`sudo`) rights using Conda, Mamba, or Pixi. All compilers, libraries, and Python modules are installed inside your home directory.

---

## 1. Setting Up the User-Space Environment

### Step 1: Install Miniforge (if you do not have Conda/Mamba)
```bash
# Download and install Miniforge into ~/miniforge3 (zero root required)
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-$(uname -m).sh"
bash Miniforge3-Linux-$(uname -m).sh -b -p $HOME/miniforge3
source $HOME/miniforge3/bin/activate
```

*(Alternatively, if using [Pixi](https://pixi.sh): `curl -fsSL https://pixi.sh/install.sh | bash`).*

### Step 2: Create an Isolated Environment with C++23 Compiler & Libraries
Conda-forge provides modern C++23 compilers (Clang 19 with libc++) and all required math and data libraries:

```bash
conda create -n jgap-env -c conda-forge -y \
    clangxx=19 llvm-openmp \
    cmake ninja \
    hdf5 openblas \
    python=3.12 pybind11 numpy scipy ase

conda activate jgap-env
```

---

## 2. Configure, Build, and Install

Point CMake to the compiler and installation prefix within the active Conda environment (`$CONDA_PREFIX`):

```bash
# 1. Configure Release build targeting the active conda environment
cmake --preset release \
      -DCMAKE_C_COMPILER=$(which clang) \
      -DCMAKE_CXX_COMPILER=$(which clang++) \
      -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX \
      -DPython3_EXECUTABLE=$(which python)

# 2. Build and install library, CLI tools, and Python bindings
cmake --build --preset install
```

---

## 3. Verification

### Verify Python & ASE Integration
```bash
python -c "import jgap; print('JGAP loaded successfully:', jgap.__file__)"
```

### Verify Standalone C++ Compilation (Without CMake)
Because all headers and libraries are located in `$CONDA_PREFIX`, use the active compiler:
```bash
clang++ -std=c++23 -O3 -march=native \
    -I$CONDA_PREFIX/include \
    -L$CONDA_PREFIX/lib \
    -Wl,-rpath,$CONDA_PREFIX/lib \
    examples/basic_fit/BasicFit.cpp -ljgap -o basic_fit

./basic_fit test/resources/structure-databases/db_Fe.xyz fe_test --ram-limit 2.0
```
