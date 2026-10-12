# Installation Guide: Linux (No `sudo` — Conda / Mamba / Pixi)

This guide covers setting up `jgap` dependencies entirely in user space without administrator (`sudo`) rights using Conda.

---

## 1. Setting Up the User-Space Environment

### Step 1: Install Conda (if you do not have it yet)
If you do not have Conda installed, you can find the installation instructions [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).

### Step 2: Create and Activate Environment
```bash
conda init # if it has never been initialized before, then reload the terminal or run `source ~/.bashrc`  
conda create -n jgap-env -y
conda activate jgap-env
```

### Step 3: Install Dependencies (Select Only What You Need)

You can copy and paste only the components relevant to your system:

#### A. C++23 Compiler (GCC & G++)

Check that your C++ compiler supports C++23. 
You can do this by running, e.g. `g++ --version`, and checking the version number. 
If it is less than 15 for G++ or 19 for Clang, you will need to install a newer version:

```bash
conda install -c conda-forge -y gcc gxx
```

#### B. Build System (CMake; Ninja is optional)
If `cmake --version` shows anything below 3.25, you will need to install a newer version:
```bash
conda install -c conda-forge -y cmake
```
If the Ninja is not installed, run:
```bash
conda install -c conda-forge -y ninja
```

#### C. Dependencies (HDF5 & OpenBLAS)
```bash
conda install -c conda-forge -y hdf5 openblas
```

#### D. Python Environment & Scientific Stack (ASE)
To avoid potential conflicts, it is recommended to install Python within conda's virtual environment:
```bash
conda install -c conda-forge -y python=3.14 pybind11 numpy scipy ase elastic
```
