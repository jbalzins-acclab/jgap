# JGAP Conventions & Guidelines

This document establishes the physical, mathematical, and coding conventions adopted throughout `jgap`.

---

## 1. Physical Units & Coordinate Systems

All internal computations and user interfaces strictly adhere to standard atomic simulation units:

| Quantity | Unit | Notes / Conventions |
| :--- | :--- | :--- |
| **Length / Positions** | Ångströms ($\text{Å}$) | Cartesian coordinates $(x, y, z)$ |
| **Energy** | Electron-volts ($\text{eV}$) | Total energy or per-atom energy |
| **Forces** | $\text{eV} / \text{Å}$ | Negative gradient of potential energy: $\mathbf{F}_i = -\nabla_{\mathbf{r}_i} E$ |
| **Virials** | $\text{eV}$ | Internal virial tensor $\mathbf{\Xi} = \sum_i \mathbf{r}_i \otimes \mathbf{F}_i$ |
| **Stress** | $\text{GPa}$ or $\text{eV} / \text{Å}^3$ | Cauchy stress tensor $\boldsymbol{\sigma} = -\frac{1}{V} \mathbf{\Xi}$ |
| **Stress Voigt Order** | — | Standard ASE ordering: $[\sigma_{xx}, \sigma_{yy}, \sigma_{zz}, \sigma_{yz}, \sigma_{xz}, \sigma_{xy}]$ |
| **Cell / Lattice** | $\text{Å}$ | Row-major $3 \times 3$ matrix representing basis vectors $(\mathbf{a}, \mathbf{b}, \mathbf{c})$ |

---

## 2. Memory Layouts & Data Representation

* **Atomic Positions & Forces**:
  Stored in `Atoms` as contiguous `std::vector<Vector3>`. This ensures vectorized cache alignment and zero-copy conversion to and from NumPy arrays of shape `(N, 3)` with `dtype=float64`.
* **Per-Atom Extra Quantities**:
  Stored in `Atoms::extra_arrays` using `Matrix<RowMajor, T>`. Supports dynamic row appending (`appendRow`) and removal (`removeRow`) when atoms are dynamically added or removed.
* **Lattice Matrices**:
  $3 \times 3$ matrices stored in row-major order:
  $$\mathbf{R} = \begin{pmatrix} a_x & a_y & a_z \\ b_x & b_y & b_z \\ c_x & c_y & c_z \end{pmatrix}$$
* **Linear Algebra Matrices**:
  Core linear algebra matrices interfacing with Eigen and BLAS default to `Eigen::MatrixXd` (ColMajor) for BLAS compatibility or `Matrix<RowMajor, double>` where row-streaming is required.

---

## 3. C++ Coding & Design Standards

* **C++23 Standard**:
  Code leverages modern C++23 features (concepts, designated initializers, `std::span`, and `#embed` where supported).
* **Value Semantics with `ValuePtr<T>`**:
  Polymorphic components (cutoffs, kernels, transformations) are managed via `ValuePtr<T>`, which provides deep-copying value semantics via a virtual `.clone()` method while preventing pointer slicing.
* **Code Formatting**:
  All C++ files follow the rules specified in the root [`.clang-format`](file:///Users/jegorsbalzins/jgap/.clang-format). Run `clang-format -i <file>` prior to submitting changes.
* **Shared Library Symbol Visibility**:
  `jgap_lib` is strictly built as a shared library with default symbol visibility to support runtime constructor-based registration (`REGISTER_SERIALIZATION`).

---

## 4. Concurrency & Thread-Safety

* **OpenMP Parallelization (`HAS_OPENMP`)**:
  Loops over atom pairs, neighbor lists, configurations, and grid evaluation use OpenMP (`#pragma omp parallel for`).
* **Default Thread Count**:
  By default, OpenMP automatically executes with the maximum available logical cores on the system (`omp_get_num_procs()`).
* **Environment Variable**:
  Users can restrict or throttle CPU thread allocation by setting `OMP_NUM_THREADS` (e.g. to avoid CPU oversubscription with BLAS threads or when sharing compute nodes).
* **Const Correctness & Thread Safety**:
  Evaluation methods (`Potential::calculateEnergy`, `Kernel::calculateCovariance`, `CutoffFunction::evaluate`) are strictly `const` and thread-safe for concurrent read access across multiple threads or configurations.
* **Python Global Interpreter Lock (GIL)**:
  Computationally heavy C++ routines (`standard_gap_fit`, `standard_tabulation`, energy/force evaluation) explicitly release the Python GIL (`py::call_guard<py::gil_scoped_release>()`), enabling true parallel multi-threading within Python scripts.
