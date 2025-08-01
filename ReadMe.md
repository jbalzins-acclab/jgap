# $\vec{ȷ}GAP$
## Overview
### Fit 2b+3b+EAM GAP fit
- On small databases fit coefficients exactly match QUIP output(without virial fit; see /test)
- Slightly different "uniform" sparsification is available
- Screened coulomb pre-fit for all element pairs(see resources/dmol-screening-fit & */core/potentials/ZblPotential.cpp)
- Significant speedup in kernel matrix formation & RAM usage improvement compared to QUIP with basic compilation:
  - ~20 sec(on my laptop) for Iron potential | ~500Mb
  - ~20 sec(on Puhti node; ~15 min total(+QR..) on a desktop but tested a while ago) for FeNi potential | ~7Gb RAM
  - 3 min(on Puhti node) for CrMnFeNi potential | ~110Gb(shown with "seff", but allocation failed when 120Gb were reserved on last fit attempt) with virial fit
  - more to be tested.
  - RAM usage can be estimated from logs(look for matrix size).
  - ! Linear algebra is slower than QUIP for now (2.2h => 4h for CrMnFeNi)
- Output in the json format not compatible with QUIP - separate app is compiled to use it.
- Per config-type regularization not implemented yet, but $\sigma$'s can be specified per structure in ext-xyz(see Utils.cpp)
### Tabulate 2b+3b+EAM 
- Very fast(around a minute on my laptop to tabulate CrMnFeNi)
- Output in .tabgap+.eam.fs (I'm not sure if non-2b+3b+EAM works correctly)

## Compilation/Run guide
### Prerequisites 
- CMake 3.26+
- c++ complier supporting c++23(c++20 might be enough, but that wasn't tested in a while)
  - do "module load gcc/14.2" on Puhti
- VCPKG:
  - (Follow instructions: https://learn.microsoft.com/en-gb/vcpkg/get_started/get-started?pivots=shell-bash))
    - git clone https://github.com/microsoft/vcpkg.git
    - cd vcpkg && ./bootstrap-vcpkg.sh
    - nano ~/.bashrc or ~/.zshrc
    - export VCPKG_ROOT=/path/to/vcpkg
    - export PATH=$VCPKG_ROOT:$PATH
  - in root project dir: vcpkg install
  - (vcpkg integrate to see what to add to cmake params)
### Compile
- Run something like: 
  - cmake -B build -DCMAKE_TOOLCHAIN_FILE=$VCPKG_ROOT/scripts/buildsystems/vcpkg.cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-O3 -DNDEBUG"
  - on Puhti: cmake -B build-test   -DCMAKE_TOOLCHAIN_FILE=/users/balzinsj/vcpkg/scripts/buildsystems/vcpkg.cmake   -DCMAKE_CXX_COMPILER=g++  -DCMAKE_CXX_FLAGS="-g -O3 -march=native -ffast-math -funroll-loops -mprefer-vector-width=512" 
- cmake --build build -j ...
- This should produce 3 executables:
  - jgap_fit_app - GAP fitting
  - jgap_predict_app - to use the GAP potential
  - jgap_tabulate_app - tabulate GAP potential
### Run
(see param samples in /resources)
- jgap_fit_app fit_param_file.json => outputs potential.json
- jgap_predict_app potential.json input.xyz output.xyz
- jgap_tabulate_app tabulation_params.json