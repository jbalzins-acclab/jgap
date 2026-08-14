# JGAP
## Overview


## Compilation/Run guide
### Prerequisites
- CMake 3.11+
- A C++23 compiler. **GCC 15+ is recommended**: jgap uses `#embed` (GCC 15 / Clang ≥ 19) to bake the
  built-in screening datasets into the binary. Older C++23 compilers (GCC 14, AppleClang) still build
  and fit — without `#embed` those datasets are instead read at runtime from
  `resources/dmol-screening-fit/`, so you must run with the `resources/` folder present (or pass an
  explicit screening dataset file; see `standard_fit`'s optional argument and `StandardGapParams::screened_coulomb_dataset_file`).
  Get a suitable compiler if needed:
  - Puhti / HPC modules: `module load gcc/15` (or the newest available)
  - Homebrew: `brew install gcc` then `export CC=gcc-15 CXX=g++-15` (the default `c++` on macOS is
    AppleClang, so you must point CMake at the Homebrew GCC — via these vars or
    `-DCMAKE_CXX_COMPILER=g++-15`)
  - conda (no sudo): `conda install -c conda-forge 'gxx>=15'` then
    `export CC=$CONDA_PREFIX/bin/gcc CXX=$CONDA_PREFIX/bin/g++`
  - Spack: `spack install gcc@15` (then `spack load gcc@15`)

### Dependencies
- **Required:** Eigen3 (header-only), oneTBB, HighFive (+ HDF5).
- **Optional:** pugixml (QUIP `.xml` conversion — see below); a BLAS such as OpenBLAS.
  - A BLAS is **strongly recommended** — it does the heavy linear algebra and is a large speedup for
    fitting. Without one, jgap falls back to Eigen's own routines and enables Eigen's OpenMP
    multi-threading (if OpenMP is available) so they at least run in parallel; a BLAS is still faster.

Install them with **any one** of the package managers below, then point CMake at them:
- vcpkg → pass its toolchain file (`-DCMAKE_TOOLCHAIN_FILE=$VCPKG_ROOT/scripts/buildsystems/vcpkg.cmake`)
- everything else → pass the install prefix(es) via `-DCMAKE_PREFIX_PATH=...`

CMake finds each dependency through its CMake config package, so any source that installs those works.
A missing optional dependency is detected automatically (a `-- ...` status line reports it).

#### vcpkg (manifest in `vcpkg.json`)
> Recommended on personal machines (laptops/desktops), where it gives a clean, reproducible build.
> On HPC clusters it works but use it with care: vcpkg builds every dependency from source and keeps
> large build trees and caches, creating a great many files — which can strain quota- or inode-limited
> shared filesystems. There, prefer the cluster's environment modules, Spack, or conda.

```bash
git clone https://github.com/microsoft/vcpkg.git
cd vcpkg && ./bootstrap-vcpkg.sh
export VCPKG_ROOT=$PWD && export PATH=$VCPKG_ROOT:$PATH
cd <jgap>   # vcpkg install reads vcpkg.json
vcpkg install
# then configure jgap with -DCMAKE_TOOLCHAIN_FILE=$VCPKG_ROOT/scripts/buildsystems/vcpkg.cmake
```

#### Homebrew (macOS / Linux)
```bash
brew install cmake eigen tbb hdf5 openblas pugixml
# HighFive is not in homebrew-core — install it via conda/spack or from git (below).
# configure with: -DCMAKE_PREFIX_PATH="$(brew --prefix)"
```

#### conda (conda-forge)
```bash
conda install -c conda-forge eigen tbb-devel highfive hdf5 openblas pugixml
# configure with: -DCMAKE_PREFIX_PATH="$CONDA_PREFIX"
```

#### Spack
```bash
spack install eigen intel-tbb highfive openblas pugixml
spack load eigen intel-tbb highfive openblas pugixml   # sets CMAKE_PREFIX_PATH
```

#### From git (build/install from source)
Useful for header-only or unpackaged deps. Each ships a CMake install; install to a common `<prefix>`
and add it to `CMAKE_PREFIX_PATH`. Example (HighFive, which needs HDF5):
```bash
git clone --depth 1 https://github.com/BlueBrain/HighFive
cmake -S HighFive -B HighFive/build -DHIGHFIVE_UNIT_TESTS=OFF -DHIGHFIVE_EXAMPLES=OFF \
      -DCMAKE_INSTALL_PREFIX=<prefix>
cmake --build HighFive/build --target install
```
Eigen (https://gitlab.com/libeigen/eigen) and pugixml (https://github.com/zeux/pugixml) install the same
way; oneTBB is at https://github.com/uxlfoundation/oneTBB.

### Optional: XML (QUIP) support
QUIP `.xml` conversion needs `pugixml`. If `pugixml` is not found, CMake automatically drops the
`QuipXmlConverter` sources and the `jgap_convert` app (a status line reports which way it went);
everything else builds unchanged.
- Other package managers: simply install (or don't install) `pugixml`.
- vcpkg: `pugixml` is the default-on `xml` feature; build without it via
  `vcpkg install --no-default-features` or, at configure time,
  `-DVCPKG_MANIFEST_NO_DEFAULT_FEATURES=ON`.
### Compile
Configure with whichever dependency source you used (see Dependencies), then build:
```bash
# with vcpkg:
cmake -B build -DCMAKE_TOOLCHAIN_FILE=$VCPKG_ROOT/scripts/buildsystems/vcpkg.cmake \
      -DCMAKE_BUILD_TYPE=Release

# without vcpkg (Homebrew / conda / Spack / from-source) — point at the install prefix(es) instead:
cmake -B build -DCMAKE_PREFIX_PATH="$(brew --prefix)" -DCMAKE_BUILD_TYPE=Release   # e.g. Homebrew
# (conda: $CONDA_PREFIX; Spack: usually set by `spack load`; multiple prefixes: ";"-separated)

cmake --build build -j
```
Add `-DCMAKE_CXX_COMPILER=g++-15` (or your C++23 compiler) if it isn't the default.
### Install (to use jgap as a library)
To consume jgap from another CMake project (see `examples/`), install it to a prefix:
```bash
cmake --install build --prefix <prefix>
```
A downstream project then only needs `find_package(jgap CONFIG REQUIRED)` and
`target_link_libraries(<target> PRIVATE jgap::jgap)` — that pulls in jgap's headers and its public
dependencies. (The dependency CMake configs must be on `CMAKE_PREFIX_PATH`; when jgap is built with
vcpkg they live in `<build>/vcpkg_installed/<triplet>`.)
This produces:
- `libjgap` — the library (link it as `jgap::jgap`; see "Install" above and `examples/`).
- `jgap` — CLI for using a serialized potential.
- `jgap_convert` — QUIP `.xml` → `.h5` converter (only when `pugixml` is available).
### Run
- `jgap --predict <pot.h5> <in.xyz> <out.xyz>` — predict energy/forces for every frame in `in.xyz`.
- `jgap --tabulate <pot.h5> [r_min_3b=..] [max_eam_density=..] [n_grid_2b=..] [n_grid_3b=a,b,c]` —
  writes `<pot>.tabgap.h5` plus the EAM `<pot>.eam.fs` file(s).
- `jgap_convert <quip.xml> [out.h5]` — convert a QUIP potential to a serialized `.h5`.
  - warn: sensitive to format changes in quip.xml — check the code upon error.

Fitting a potential is done through the library API (no dedicated CLI yet) — see `examples/Main.cpp`.