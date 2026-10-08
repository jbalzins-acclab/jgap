# Installation Guide: Linux (No `sudo` — Apptainer / Singularity `.sif`)

This guide covers running and building `jgap` on **any Linux distribution** (Ubuntu, RHEL, Rocky Linux, CentOS, Debian, Fedora, Arch) without administrator (`sudo`) access using **Apptainer** (formerly Singularity) or a pre-built `.sif` image.

---

## 1. Why Apptainer / `.sif`?

* **Zero `sudo` or root privileges**: Unlike Docker, Apptainer runs entirely under your standard unprivileged user account.
* **Universal Linux compatibility**: The `.sif` image bundles its own Ubuntu 24.04 environment, GCC 15, HDF5, OpenBLAS, and Python. It runs identically across any Linux kernel version and distro without glibc or library conflicts.
* **Native performance**: Apptainer does not emulate a virtual machine; it uses direct host Linux namespaces with full CPU vectorization (`-march=native`) and OpenMP multi-threading.
* **Automatic filesystem binding**: Your `$HOME` directory and current working directory (`$PWD`) are automatically mounted inside the container.

---

## 2. Obtaining or Building `jgap.sif`

### Option A: Pull Pre-built Container (Recommended)
If a pre-built image is hosted on a container registry:
```bash
# Pull or convert directly in user space (no root required)
apptainer pull jgap.sif docker://ghcr.io/jbalzins-acclab/jgap:latest
```

### Option B: Build from Definition File (`containers/jgap.def`)
Because building a new `.sif` image from a definition file requires root during the package installation step, build it on your personal computer (or via Docker / GitHub Actions) and transfer `jgap.sif` to your Linux server:

```bash
# On a machine with root/sudo or Apptainer fakeroot:
apptainer build jgap.sif containers/jgap.def

# Transfer to your remote server:
scp jgap.sif user@cluster.institution.edu:~/
```

---

## 3. Using the `.sif` Container on the Target Machine

Once `jgap.sif` is on your machine, you can use it immediately without any installation:

### A. Interactive Development Shell
Enter an interactive shell with GCC 15, CMake, Ninja, and Python all available:
```bash
apptainer shell jgap.sif
```

Inside the shell, your local repository files are accessible directly:
```bash
# 1. Configure and build jgap
cmake --preset release
cmake --build --preset install

# 2. Compile standalone examples without CMake
g++ -std=c++23 -O3 -march=native examples/basic_fit/BasicFit.cpp -ljgap -o basic_fit

# 3. Run fit
./basic_fit test/resources/structure-databases/db_Fe.xyz fe_pot --ram-limit 2.0
```

### B. Direct Command Execution (`apptainer exec`)
You do not need to open an interactive shell; you can run any tool directly through the container:

```bash
# Run C++ fit executable:
apptainer exec jgap.sif ./basic_fit test/resources/structure-databases/db_Fe.xyz fe_pot --ram-limit 2.0

# Run Python potential evaluation with ASE:
apptainer exec jgap.sif python3 examples/ase_integration/test_potential_with_ase.py \
    fe_pot.tabgap.h5 \
    test/resources/structure-databases/db_Fe.xyz
```

### C. Submitting Batch Jobs on HPC (SLURM / PBS)
Apptainer works seamlessly in HPC job scripts:

```bash
#!/bin/bash
#SBATCH --job-name=jgap_fit
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --time=02:00:00

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Execute directly through the container
apptainer exec jgap.sif ./standard_fit data/train.xyz fe_pot --ram-limit 16.0
```
