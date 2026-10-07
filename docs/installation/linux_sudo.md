# Installation Guide: Linux (with `sudo` rights)

This guide covers building and installing `jgap` natively on Linux when you have administrative (`sudo`) access to install packages.

---

## 1. Prerequisites: GCC 15+ & System Libraries

`jgap` utilizes C++23 features (including standard multi-dimensional spans and language extensions) that require **GCC $\ge$ 15** or **Clang $\ge$ 19**.

### Ubuntu / Debian (`apt`)
On Ubuntu 24.04 and earlier, install GCC 15 from the official Ubuntu Toolchain PPA:

```bash
# 1. Add Ubuntu Toolchain PPA for GCC 15
sudo apt update && sudo apt install -y software-properties-common
sudo add-apt-repository -y ppa:ubuntu-toolchain-r/test
sudo apt update

# 2. Install GCC 15, CMake, Ninja, HDF5, OpenBLAS, and Python dev packages
sudo apt install -y --no-install-recommends \
    gcc-15 g++-15 cmake ninja-build \
    libhdf5-dev libopenblas-dev \
    python3 python3-dev python3-pip python3-venv git
```

### Fedora / RHEL (`dnf`)
```bash
sudo dnf install -y gcc-c++ cmake ninja-build \
                    hdf5-devel openblas-devel \
                    python3-devel python3-pip
```

### Arch Linux (`pacman`)
```bash
sudo pacman -S gcc cmake ninja hdf5 openblas python python-pip
```

---

## 2. Setting Up Python Environment

```bash
# 1. Create and activate a virtual environment
python3 -m venv ~/jgap_env
source ~/jgap_env/bin/activate

# 2. Install pybind11 and scientific packages
pip install --upgrade pip
pip install pybind11 numpy scipy ase
```

---

## 3. Configure, Build, and Install

Pass `gcc-15` and `g++-15` explicitly to CMake, pointing `CMAKE_INSTALL_PREFIX` to your active virtual environment:

```bash
# 1. Configure optimized Release build
cmake --preset release \
      -DCMAKE_C_COMPILER=gcc-15 \
      -DCMAKE_CXX_COMPILER=g++-15 \
      -DPython3_EXECUTABLE=$(which python) \
      -DCMAKE_INSTALL_PREFIX=$VIRTUAL_ENV

# 2. Build and install library, CLI tools, and Python extension
cmake --build --preset install
```

> [!NOTE]
> If installing system-wide or to user space outside a virtual environment, omit `CMAKE_INSTALL_PREFIX` (defaults to `$HOME/.local`). If installing to `$HOME/.local`, make sure `$HOME/.local/lib` is in your `LD_LIBRARY_PATH`:
> ```bash
> export LD_LIBRARY_PATH="$HOME/.local/lib:$LD_LIBRARY_PATH"
> ```

---

## 4. Verification

### Run Unit Tests & Validation Benchmarks
```bash
# Unit tests
cmake --workflow --preset ci

# Validation suite vs reference potentials
cmake --workflow --preset validation
```

### Verify Standalone C++ Compilation (Without CMake)
```bash
g++-15 -std=c++23 -O3 -march=native examples/basic_fit/BasicFit.cpp -ljgap -o basic_fit
./basic_fit test/resources/structure-databases/db_Fe.xyz fe_test --ram-limit 2.0
```

### Verify Python & ASE Integration
```bash
python -c "import jgap; print('JGAP loaded successfully:', jgap.__file__)"
```
