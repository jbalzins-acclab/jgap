# Benchmarking & Validation Setup

[Explanation of the benchmarking and validation process goes here. This section will be filled in later to describe the workflow, speed comparisons, test cases, and goals.]

## Prerequisites

Before running the installation scripts, ensure your system has the required dependencies installed. You will need standard development tools, Fortran compilers, CMake, Git, linear algebra libraries (LAPACK/BLAS), and OpenMP support.

### macOS (Apple Silicon or Intel)
If you are on macOS, the easiest way to get the dependencies is via Homebrew:

```sh
# Install compilers, cmake, OpenMP, and LAPACK/BLAS
brew install gcc cmake libomp lapack openblas
```

### Ubuntu / Debian Linux
If you are on a Debian-based system, you can use `apt`:

```sh
sudo apt update
sudo apt install build-essential gfortran cmake git liblapack-dev libblas-dev libomp-dev
```

## Running the Benchmarks

To automatically download benchmarking datasets (gap-data) and install all dependencies (QUIP, TabGAP, and LAMMPS - took about 6 minutes to run on base MacBook Air M2)
and execute the benchmarks, simply run the main script from this directory:

```sh
chmod +x run_benchmarking.sh
./run_benchmarking.sh
```

Alternatively, you can run the individual installation scripts directly if you only need to rebuild a specific component:
* `./download_gap_data.sh`
* `./install_quip.sh`
* `./install_tabgap.sh`
* `./install_lammps.sh`
