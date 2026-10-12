# Installation Guide: macOS (Apple Silicon with Homebrew)

This guide covers installation of `jgap` dependencies on Apple Silicon macOS (M-series) using AppleClang and Homebrew.

---

## 1. Prerequisites & Dependencies

On macOS, Homebrew is assumed to be available. AppleClang (via Xcode Command Line Tools) provides full C++23 language support. Apple Accelerate is included with macOS and provides native, highly optimized BLAS acceleration (OpenBLAS is not needed).

### Step 1: Install Xcode Command Line Tools (if not already installed)
```bash
xcode-select --install
```
*(Verify by running `clang++ --version`; AppleClang $\ge$ 19 is supported).*

### Step 2: Install Libraries via Homebrew

You can install only the packages you need:

If CMake >= 3.25 is not installed:
  ```bash
  brew install cmake
  ```
If Ninja is not installed:
  ```bash
  brew install ninja
  ```

If HDF5 Database is not installed:
  ```bash
  brew install hdf5 
  ```
---

## 2. Setting Up Python Environment

To use `jgap`'s Python interface and ASE calculator, activate a Python virtual environment:

```bash
# 1. Create and activate a virtual environment
python3 -m venv ~/my_env
source ~/my_env/bin/activate
pip install --upgrade pip
```

And install Python dependencies:
```bash
  pip install pybind11 numpy scipy ase elastic
```

---