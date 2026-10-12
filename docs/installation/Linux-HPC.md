# Installation Guide: Linux on HPC Clusters

This guide covers setting up `jgap` compilation instructions on HPC clusters where `sudo` is not available. 
On many HPC clusters, a lot of libraries come pre-installed and are loaded using the `module` command.  
This is not always the case, so check if the required libraries are installed on your cluster via `module avail` or `module spider <module_name/[version]>`.

---

## 1. Loading the Environment

To load the dependencies that we assume to be preinstalled:

   ```bash
   module purge
   module load gcc/15     # Or llvm/clang-19
   module load openblas
   module load lapack     # For LAPACKE QR acceleration (if separated from OpenBLAS)
   module load hdf5
   module load cmake
   ```

### 2. (Optional) Create/Activate Python Virtual Environment

```bash
# 1. Create and activate a virtual environment
python3 -m venv ~/my_env
source ~/my_env/bin/activate
pip install --upgrade pip
```

### 3. Install Ninja if not available

Ninja is optional (load if available, or install via 'pip install ninja')
```bash
module load ninja 2>/dev/null || true
```

### 4. Install Python Dependencies if not available

```bash
pip install pybind11 numpy scipy ase elastic
```