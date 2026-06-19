# jgap examples

Small standalone programs that use the jgap library. Each is its own CMake project that just does
`find_package(jgap)` + links `jgap::jgap`.

- **`standard_fit/`** — `standard_fit <training.xyz> <output_prefix>`: fits a 2b+3b+EAM GAP, serializes it
  to `<output_prefix>.jgap.h5`, and tabulates it to `<output_prefix>.tabgap.h5` + `<output_prefix>.eam.fs`
  (logging where they were written).
- **`read_and_predict/`** — `read_and_predict <pot.h5> <in.xyz> <out.xyz>`: loads any serialized potential
  and writes per-frame predictions.
- **`CustomFit/`** — `custom_fit [training.xyz] [output_prefix]`: a hand-built FeNi fit (instead of
  `standardGapFit`) showing per-component flexibility — a different EAM pair function per element pair
  (FSGen / Coscutoff / Polycutoff) and a per-triplet 3-body energy scale that grows with the Ni count.
  Baseline numbers are taken from `test-local/fit-params.json` (inlined, as jgap has no JSON reader).
- **`new_app.sh`** — scaffolds a new app (creates `<AppName>/<AppName>.cpp` + a `CMakeLists.txt` and prints
  build/run instructions).

## Prerequisite

Build jgap first — see the top-level `ReadMe.md` (which also covers getting the dependencies via vcpkg /
Homebrew / conda / Spack / git). **No install step is needed**: jgap exports its build tree, so the
examples can find it straight from the jgap build directory.

You will point CMake at two prefixes:
- `<jgap_build>` — the jgap build directory (holds `jgapConfig.cmake`)
- `<deps_prefix>` — where jgap's dependency CMake configs live: with vcpkg this is
  `<jgap_build>/vcpkg_installed/<triplet>`; otherwise the manager prefix (`$(brew --prefix)`,
  `$CONDA_PREFIX`, a Spack/from-source prefix, ...)

The examples don't use any package manager themselves.

## Build and run an example

```bash
cd standard_fit
cmake -B build -DCMAKE_PREFIX_PATH="<jgap_build>;<deps_prefix>"
cmake --build build
./build/standard_fit <training.xyz> fe-pot       # e.g. ../../resources/xyz-samples/db_Fe.xyz

cd ../read_and_predict
cmake -B build -DCMAKE_PREFIX_PATH="<jgap_build>;<deps_prefix>"
cmake --build build
./build/read_and_predict fe-pot.jgap.h5 <in.xyz> out.xyz
```

Each example's whole `CMakeLists.txt` is just:

```cmake
find_package(jgap CONFIG REQUIRED)
add_executable(<name> <Name>.cpp)
target_link_libraries(<name> PRIVATE jgap::jgap)
```

## Scaffold a new app

```bash
./new_app.sh
```
It asks for the app name, the jgap build directory, and whether jgap was built with vcpkg (to locate the
dependency configs), then generates `<AppName>/<AppName>.cpp` and `<AppName>/CMakeLists.txt` and prints
the build/run commands.
