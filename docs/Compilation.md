# Compilation & Installation Guide

`jgap` is a high-performance C++23 library, CLI toolkit, and Python framework. Because different environments have different privileges and package managers, dedicated installation guides are provided below.

---

## 1. Quick Selector: Choose Your Environment

Select the guide matching your platform and access level:

| Platform | Access Level / Tool | Guide |
| :--- | :--- | :--- |
| **macOS** (Apple Silicon) | Homebrew available | [**macOS Installation Guide**](installation/macos.md) |
| **Linux** (Ubuntu, Debian, Fedora, Arch) | `sudo` root rights available | [**Linux (`sudo`) Guide**](installation/linux_sudo.md) |
| **Linux** (Any distribution / Server) | No `sudo` — Conda / Mamba / Pixi | [**Linux (Conda / User-Space) Guide**](installation/linux_conda.md) |
| **Linux** (Any distribution / HPC) | No `sudo` — Apptainer / Singularity `.sif` | [**Linux (Apptainer `.sif`) Guide**](installation/linux_apptainer.md) |

---

## 2. Core Technical Requirements Summary

* **Compilers**:
  * **Linux**: GCC $\ge$ 15 (required for GCC) or Clang $\ge$ 19.
  * **macOS**: AppleClang $\ge$ 16 (Xcode 16+ on Apple Silicon).
* **Build System**: CMake $\ge$ 3.25 and Ninja.
* **C++ Libraries**: HDF5 (`libhdf5`), OpenBLAS (or Apple Accelerate), OpenMP.
* **Python**: Python $\ge$ 3.8 with `pybind11`, `numpy`, and `ase`.

---

## 3. CMake Workflow Presets

`jgap` provides unified CMake presets (`CMakePresets.json`) across all platforms:

| Preset Name | Type | Description |
| :--- | :--- | :--- |
| `release` | Configure | Optimized release configuration (`-O3 -ffast-math -march=native`). |
| `release-tests` | Configure | Release build with unit tests enabled (`JGAP_BUILD_TESTS=ON`). |
| `debug` | Configure | Debug symbols, compile commands export (`compile_commands.json`). |
| `validation` | Configure | Release configuration with automated validation test suite. |
| `asan` | Configure | AddressSanitizer and UndefinedBehaviorSanitizer instrumentation. |
| `install` | Build | Builds Release targets and installs them to `CMAKE_INSTALL_PREFIX`. |

### One-Command Workflow Execution
```bash
# Run unit tests in Release mode
cmake --workflow --preset ci

# Build and install Release binaries
cmake --workflow --preset install

# Build and execute validation test suite
cmake --workflow --preset validation

# Debug build and test execution
cmake --workflow --preset dev
```

### AddressSanitizer Run
```bash
cmake --preset asan
cmake --build --preset asan
ctest --preset asan
```

---

## 4. Troubleshooting & FAQs

### Python Extension Links to Wrong Python Version
If CMake selects a system Python rather than your active virtual environment's Python:
1. Delete the `build/` directory: `rm -rf build`.
2. Explicitly specify the interpreter path during configuration:
   ```bash
   cmake --preset release -DPython3_EXECUTABLE=$(which python)
   ```

### Pybind11 Not Found During CMake Configuration
If you see `-- pybind11 not found: skipping _jgap Python module build`:
Ensure `pybind11` is installed inside your active environment:
```bash
pip install pybind11
```
and re-run CMake configuration.
