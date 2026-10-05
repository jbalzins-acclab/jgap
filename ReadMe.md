# JGAP

A high-performance C++23 library, command-line toolkit, and Python framework for fitting, tabulating, and evaluating **Gaussian Approximation Potentials (GAP)** and **Tabulated GAP (tabGAP)**.

---

## 1. Quickstart: Compile & Install

You can compile and install `jgap` (both the C++ shared library, CLI tools, and the Python interface into your active virtual environment) in just two commands:

```bash
# 1. Configure optimized Release build targeting your active Python virtual environment
cmake --preset release -DPython3_EXECUTABLE=$(which python) -DCMAKE_INSTALL_PREFIX=$VIRTUAL_ENV

# 2. Build and install library, CLI executables, and Python package
cmake --build --preset install
```

> [!NOTE]
> If installing to your user directory outside a virtual environment, use `-DCMAKE_INSTALL_PREFIX=$HOME/.local` (or run `cmake --workflow --preset install`). For Conda environments, use `-DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX`.
>
> For full prerequisites (C++23 compilers, HDF5, OpenBLAS), package manager commands (macOS, Ubuntu, Conda), HPC cluster instructions, and troubleshooting, see the **[Compilation & Installation Guide](docs/Compilation.md)**.

---

## 2. Python Quickstart & ASE Calculator

`jgap` integrates directly with Python and the [Atomic Simulation Environment (ASE)](https://wiki.fysik.dtu.dk/ase/):

```python
import jgap
from jgap.ase import JGAPCalculator
from ase.io import read

# 1. Load training database
training_data = jgap.read_atoms("train.xyz")

# 2. Fit a standard 2-body + 3-body + EAM GAP potential
params = jgap.StandardGapParams(
    seed=120,
    n_sparse3=500,
    eam_mode=jgap.EamMode.Blind
)
sigmas = jgap.PerConfigTypeRegularizationRules(
    jgap.PerConfigTypeSigmas(0.001, 0.05, 0.1, 0.02)
).determine_for_all(training_data)

jgap.standard_gap_fit("potential.jgap.h5", training_data, sigmas, params)

# 3. Tabulate into fast 3D cubic B-spline table and EAM (.tabgap.h5 and .eam.fs)
jgap.standard_tabulation("potential.jgap.h5", "potential")

# 4. Evaluate using the ASE Calculator
atoms = read("structure.xyz")
atoms.calc = JGAPCalculator("potential.jgap.h5") # Supports .jgap.h5 and .tabgap.h5

energy = atoms.get_potential_energy()
forces = atoms.get_forces()
stress = atoms.get_stress()
```

---

## 3. Command-Line Tools

* **Predict Energy & Forces**:
  ```bash
  jgap --predict potential.jgap.h5 input.xyz output.xyz
  ```
* **Tabulate Potential (to EAM `.eam.fs` and tabGAP `.tabgap.h5`)**:
  ```bash
  jgap --tabulate potential.jgap.h5
  ```
* **Convert Legacy QUIP XML to HDF5**:
  ```bash
  jgap_convert potential.xml potential.jgap.h5
  ```

---

## 4. Multi-Threading & Performance

`jgap` utilizes OpenMP multi-threading across neighbor-list construction, energy/force evaluation, and spline tabulation.

By default, OpenMP automatically runs with the **maximum number of available CPU cores/threads** on your system. If you need to restrict or throttle core allocation (e.g. on shared HPC nodes, in CI, or when running multiple jobs in parallel), set the standard environment variable:

```bash
export OMP_NUM_THREADS=8
```

### Performance Benchmarks & Scaling
> [!NOTE]
> *(Placeholder: Detailed performance comparison plots against QUIP, evaluation latency scaling across atom counts, and memory footprint comparisons will be published here).*

---

## 5. LAMMPS Integration

* **EAM Potentials**: The generated `.eam.fs` files can be directly evaluated inside LAMMPS via standard `pair_style eam/fs`.
* **Tabulated GAP Potentials (`.tabgap.h5`)**: Direct evaluation of multi-body B-spline `.tabgap.h5` potentials inside LAMMPS simulations is supported via the external tabGAP LAMMPS package available at **[gitlab.com/jezper/tabgap](https://gitlab.com/jezper/tabgap)**.

---

## 6. Documentation Index

* **[Compilation & Installation Guide](docs/Compilation.md)**: Complete guide on toolchains, dependencies, virtual environments, HPC modules, and CMake presets.
* **[Developer & Extension Guide](docs/Development.md)**: Architectural diagrams, HDF5 serialization mechanisms, adding custom extensions (cutoffs, kernels, potentials), and running tests.
* **[Conventions & Standards](docs/Conventions.md)**: Physical units, memory layout, Voigt stress conventions, and coding style.
* **[Examples Catalog](examples/ReadMe.md)**: Guide to C++ and Python examples, including standalone single-file C++ compilation commands.

---

## 7. References & Citations

> [!NOTE]
> *(Placeholder: Bibliographic references and citations for Gaussian Approximation Potentials, tabGAP, and foundational methods will be listed here).*
