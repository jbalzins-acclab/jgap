# Compilation & Installation Guide

`jgap` is a high-performance C++23 library with a CLI toolkit and a Python interface. Since different environments have different privileges, package managers, and hardware capabilities, dedicated installation guides and workflow presets are detailed below.

---

## 1. Core Technical Requirements Summary

* **Compilers**:
  * **GCC $\ge$ 15**
  * **Clang $\ge$ 19**
* **Build System**: CMake $\ge$ 3.25 (Ninja is optional but recommended for faster builds).
* **C++ Libraries**: HDF5 (`libhdf5`), (optionally) OpenBLAS or Apple Accelerate, (optionally) OpenMP or Apple's GCD.
* **Python**: Python $\ge$ 3.9 with `pybind11`, `numpy`, and `ase` (some examples require `scipy` and `elastic`).


Some tips for installing the dependencies based on your system can be found in:

| Platform | Access Level / Tool | Guide |
| :--- | :--- | :--- |
| **macOS** (Apple Silicon) | Homebrew available | [**macOS Installation Guide**](installation/MacOS.md) |
| **Linux** | No `sudo` — Conda | [**Linux (Conda / User-Space) Guide**](installation/Linux-Conda.md) |
| **Linux on HPC Clusters** | No `sudo` - `module` command is available| [**Linux on HPC Clusters Guide**](installation/Linux-HPC.md) |

Regardless of installation type, its strongly recomnended that the following is performed in a virtual environment (e.g. `python -m venv ...`) to avoid all kinds of package conflicts.

---

## 2. Basic Compilation and Optimized `pip` Installation

Regardles of wether you want to use only the basic functionality availaible in JGAP's Python interface or you want to use JGAP as a library, 
peak performance requires system-specific compilation flags, which are disabled when using the defualt `pip install jgap` command from PyPI.

Assuming you've loaded all the needed dependencies as described in the `docs/installation` folder, and activated the virtual environmen of your choice, 
JGAP can be compiled (including Python installation to your venv) via a simple command:

```bash
   cmake --workflow --preset install
```

This will install JGAP with all its functionalities to your virtual environment, and when it is active you should be able to compile standalone tools with simple commands, like: 
```bash
c++ -std=c++23 -O3 -march=native BasicFit.cpp -ljgap -o basic_fit
```

> [!NOTE]
> If you run `cmake --workflow --preset install` without using a virtual environment, the `libjgap` library will be installed into `$HOME/.local`.
> For one-line compilation commands to work, you may need to add it to your `$PATH` in `~/.zshrc` or `~/.bashrc`:
> ```bash
> # 1. CLI tools (jgap, jgap_convert)
> export PATH="$HOME/.local/bin:$PATH"
> # 2. C++ Headers (#include <jgap/...>)
> export CPLUS_INCLUDE_PATH="$HOME/.local/include${CPLUS_INCLUDE_PATH:+:$CPLUS_INCLUDE_PATH}"
> # 3. Linker search path for -ljgap at compile-time
> export LIBRARY_PATH="$HOME/.local/lib${LIBRARY_PATH:+:$LIBRARY_PATH}"
> # 4. Dynamic loader runtime search path for libjgap.dylib
> export DYLD_LIBRARY_PATH="$HOME/.local/lib${DYLD_LIBRARY_PATH:+:$DYLD_LIBRARY_PATH}"
> export DYLD_FALLBACK_LIBRARY_PATH="$HOME/.local/lib${DYLD_FALLBACK_LIBRARY_PATH:+:$DYLD_FALLBACK_LIBRARY_PATH}"
> ```
> Additionally, adding `-Wl,-rpath,$HOME/.local/lib` to the compilation command may help.

> [!TIP]
> To clean up `$HOME/.local` run:
> ```bash
> rm -rf ~/.local/include/jgap \
>        ~/.local/lib/cmake/jgap \
>        ~/.local/lib/libjgap.* \
>        ~/.local/bin/jgap* \
>        ~/.local/lib/python*/site-packages/jgap*
> ```

You can check Python package is installed correctly by running:

```bash
python -c "import jgap; print('OK')"
```

If you only want to use JGAP's Python interface, the following command might be sufficient, given all the dependancies are available as normal:

```bash
pip install jgap --no-binary jgap \
  -Ccmake.define.CMAKE_CXX_FLAGS_RELEASE="-O3 -march=native -ffast-math"
```

If you do not need Python interface, on the other hand, you may compile via:
```bash
   cmake --preset release -DJGAP_BUILD_PYTHON=OFF
cmake --build --preset install
```

---

## 5. CMake Workflow Presets Reference

If you want to modify JGAP, you may find the following CMake presets (`CMakePresets.json`) useful:

| Preset Name | Type | Description |
| :--- | :--- | :--- |
| `release` | Configure | Optimized release configuration (`-O3 -ffast-math -march=native`). |
| `release-tests` | Configure | Release build with unit tests enabled (`JGAP_BUILD_TESTS=ON`). |
| `debug` | Configure | Debug symbols, compile commands export (`compile_commands.json`). |
| `validation` | Configure | Release configuration with automated validation test suite. |
| `asan` | Configure | AddressSanitizer and UndefinedBehaviorSanitizer instrumentation. |
| `install` | Build | Builds Release targets and installs them to `CMAKE_INSTALL_PREFIX`. |

### Execution Commands
```bash
# Build and install Release binaries in a single step (Recommended)
cmake --workflow --preset install

# Run unit tests in Release mode
cmake --workflow --preset ci

# Build and execute validation test suite
cmake --workflow --preset validation

# Debug build and test execution
cmake --workflow --preset dev

# AddressSanitizer run
cmake --preset asan && cmake --build --preset asan && ctest --preset asan
```

---

## 6. Troubleshooting & FAQs

### Ninja Build Tool is Optional
Ninja is optional. If Ninja is not available on your cluster module tree or package manager:
1. You can install it into your virtual environment in seconds via `pip install ninja`.
2. Or build using standard Makefiles:
   ```bash
   cmake -B build -DCMAKE_BUILD_TYPE=Release
   cmake --build build -j$(nproc)
   ```

### Python Extension Links to Wrong Python Version
If CMake selects a system Python rather than your active virtual environment's Python:
1. Delete the `build/` directory: `rm -rf build` .
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
