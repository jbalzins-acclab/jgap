# Installation Guide: macOS (Apple Silicon with Homebrew)

This guide covers building and installing `jgap` natively on Apple Silicon macOS (M-series) using AppleClang and Homebrew for dependencies.

---

## 1. Prerequisites & Dependencies

On macOS, Homebrew is assumed to be available. AppleClang (via Xcode Command Line Tools) provides full C++23 language support. Apple Accelerate is included with macOS and provides native, highly optimized BLAS acceleration (OpenBLAS is not needed).

### Step 1: Install Xcode Command Line Tools (if not already installed)
```bash
xcode-select --install
```
*(Verify by running `clang++ --version`; AppleClang $\ge$ 16 is supported).*

### Step 2: Install Libraries via Homebrew
Homebrew runs without `sudo`:
```bash
brew install cmake ninja hdf5 libomp
```

> [!NOTE]
> Homebrew packages are installed in `/opt/homebrew`. CMake automatically detects `/opt/homebrew/include` and `/opt/homebrew/lib`. Apple Accelerate BLAS is automatically detected and used natively out-of-the-box.

---

## 2. Setting Up Python Environment

To use `jgap`'s Python interface and ASE calculator, activate your Python virtual environment and install dependencies:

```bash
# 1. Create and activate a virtual environment
python3 -m venv ~/my_env
source ~/my_env/bin/activate

# 2. Install required Python packages
pip install --upgrade pip
pip install pybind11 numpy scipy ase
```

---

## 3. Configure, Build, and Install

Target your active virtual environment by setting `CMAKE_INSTALL_PREFIX=$VIRTUAL_ENV`:

```bash
# 1. Configure optimized Release build
cmake --preset release \
      -DPython3_EXECUTABLE=$(which python) \
      -DCMAKE_INSTALL_PREFIX=$VIRTUAL_ENV

# 2. Build and install library, CLI executables, and Python package
cmake --build --preset install
```

> [!TIP]
> If installing to your user directory outside a virtual environment, use `-DCMAKE_INSTALL_PREFIX=$HOME/.local` instead.

---

## 4. Verification

### Verify Python & ASE Integration
```bash
python -c "import jgap; print('JGAP loaded from:', jgap.__file__)"
```

### Verify Standalone C++ Compilation (Without CMake)
Once installed, standalone executables link directly against `libjgap`:
```bash
c++ -std=c++23 -O3 -march=native examples/basic_fit/BasicFit.cpp -ljgap -o basic_fit
./basic_fit test/resources/structure-databases/db_Fe.xyz fe_test --ram-limit 2.0
```
