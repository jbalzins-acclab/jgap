# Developer & Extension Guide

This guide covers the internal architecture of JGAP, how to run unit tests and validation suites, how the serialization system works, and how to extend the framework with new components.

---

## 1. System Architecture & Components

JGAP models GAP interactions through a modular, component-based potential architecture. 
A `GapPotential` (Gaussian Approximation Potential) and its tabulated counterpart `TabGapPotential`—are built from an assembly of independent **Components** alongside an optional external baseline potential.

### How Components Work

Components (`GapComponent` and `TabGapComponent`) encapsulate the physical interactions, geometric representations, and mathematical kernels for a specific interaction order (such as 2-body pair distances, 3-body angular triplets, and many-body/EAM electron densities):

1. **Neighbor Extraction & Coordinate Transformations**:
   Neighbor lists (`NeighbourLists`) are constructed for atomic configurations under periodic or free boundary conditions and expanded into 2-body pairs or 3-body triplets (`Cluster2Expansion`, `Cluster3Expansion`). Concrete coordinate transformations (e.g. `PairDistanceTransformation`, `Angle3bTransformation`, `FSGenPairFunction`) map these spatial clusters into invariant descriptor coordinates.

2. **Many-Body Aggregation (`NBodyAggregator`)**:
   For many-body and EAM components (`ManyBodyGapComponent`, `EamTGComponent`), local atomic environments cannot be described by isolated pair or triplet terms alone. Instead, an `NBodyAggregator` (such as `TwoBodySum` or `ThreeBodySum`) sums neighbor cluster contributions around each central atom into an atomic descriptor vector (e.g. local electron density $\rho_i = \sum_{j} \phi(r_{ij})$), along with derivative tensors needed for force and virial evaluations.

3. **Sparsification & Representative Basis (`Sparsifier`)**:
   Rather than treating all atomic environments in the training database as Gaussian process basis centers, each component maintains a representative subset of **Sparse Points** . All descriptors present across training configurations (`getAllDescriptors`) are passed to a `Sparsifier` (such as `HistogramUniformSparsifier` or `FromFileSparsifier`), which selects an evenly distributed basis set across descriptor space (e.g., by partitioning descriptor space into multidimensional histogram bins to ensure uniform coverage and prevent clustering in over-represented equilibrium geometries).

4. **Covariance Kernels & Covariate Evaluation**:
   Covariance kernels (e.g. `SquaredExpKernel`) are independent mathematical functions. A component combines its coordinate representation with a kernel and its sparse basis  to compute covariance blocks:
   * **$K_{MM}$**: Covariance matrix between sparse points.
   * **Covariates ($K_{NM}$)**: Evaluated per training configuration via `GapComponent::covariate(neighbour_list)`, which maps the configuration's atomic clusters (or aggregated many-body descriptors) against the sparse points to compute energy covariances and analytical force/virial derivatives for the design matrix.

5. **Linear Expansion Coefficients**:
   During fitting (`GapFit`), the linear solver solves the regularized system against DFT reference residuals ($y_{\text{target}} = y_{\text{DFT}} - y_{\text{ext}}$). Each component stores its solved linear coefficients ($c_j$) via `setCoefficients`, which scale the basis functions.

6. **Energy, Force & Virial Prediction**:
   At evaluation time, the total energy of a `GapPotential` is the sum of energies predicted across its constituent components plus any external baseline potential. Analytical derivatives with respect to atomic coordinates yield forces and virials.

7. **Tabulation Support**:
   To tabulate a GAP, low-dimensional external potentials and components evaluate their energy contributions onto discrete grids (`tabulate` / `fillTables`).

### Architectural Key Concepts
1. **Cluster Expansions & Coordinate Transformations**: Neighbor lists are expanded into 2-body pairs and 3-body triplets (`Cluster2Expansion`, `Cluster3Expansion`) before applying non-linear coordinate mappings (e.g. `Angle3bTransformation`, `FSGenPairFunction`).
2. **Many-Body Density Aggregators (`NBodyAggregator`)**: In many-body and EAM models, `TwoBodySum` or `ThreeBodySum` accumulate cluster transformations around each central atom to form atomic descriptors (such as local electron density $\rho_i$) and coordinate derivative tensors.
3. **Descriptor Sparsification (`Sparsifier`)**: Training configurations yield a dense set of descriptors. `HistogramUniformSparsifier` bins the descriptor space to select representative **Sparse Points**  with uniform coverage, preventing over-representation of common configurations.
4. **Independent Kernels & Component Assembly**: Covariance kernels are purely mathematical objects defined independently of descriptors. A GAP Component combines a specific transformation (or aggregator), its sparse basis set , and a covariance kernel.
5. **`GapPotential` as Both Input and Output of Fitting**:
   * **As Input**: `GapPotential` defines the model architecture (the collection of components, kernels, sparse points) and optionally an `optional_external_potential` (such as ZBL `ScreenedCoulombPotential` and `IsolatedAtomPotential`).
   * **Residual Fitting**: During `GapFit`, the external potential evaluates baseline contributions, which are subtracted from DFT reference targets ($y_{\text{target}} = y_{\text{DFT}} - y_{\text{ext}}$).
   * **As Output**: The linear solver (`QRGapFit` or low-memory streaming variants) solves the regularized system and assigns the resulting coefficients ($c_j$ / $\boldsymbol{\alpha}$) directly into the components of the `GapPotential`.
6. **External Simulation Engines**:
   * `.eam.fs` files run directly in LAMMPS via standard `pair_style eam/fs`.
   * Multi-body `.tabgap.h5` tables run inside LAMMPS via the external plugin repository [gitlab.com/jezper/tabgap](https://gitlab.com/jezper/tabgap).

---

## 2. Codebase Organization

The source code in `src/jgap/` maintains a clean separation between **core abstractions** and **concrete implementations**:

* **`core/`**: Foundational mathematics, interfaces, and physics data primitives:
  * `atomic/`: `Atoms`, `Lattice`, `Species`, `Virials`, neighbor lists, cluster expansions (`Cluster2Expansion`, `Cluster3Expansion`).
  * `cutoff/`: Cutoff abstract base class (`CutoffFunction`).
  * `kernels/`: Covariance kernel base interface (`Kernel`).
  * `potentials/`: Base `Potential` class, `GapPotential`, `TabGapPotential`, `CompositePotential`.
  * `sparsification/`: Sparsification base interface (`Sparsifier`).
  * `splines/`: Spline interpolation interfaces (`Spline`, `Grid`).
  * `tabulation/`: Multi-grid tabulation data structures (`TabulationData`, `TabulationParams`).
  * `transform/`: Coordinate transform interfaces (`TwoBodyTransformation`, `ThreeBodyTransformation`, `EamPairFunction`, `NBodyAggregator`).
  * `fit/`: Abstract linear solver interface (`GapFit`).
* **`impl/`**: Concrete algorithms and component implementations:
  * `cutoff/`: `CosCutoff`, `PerriotPolynomialCutoff`, `WendlandFunction`.
  * `kernels/`: `SquaredExpKernel`, `CauchyKernel`, `WendlandKernel`.
  * `fit/gap/`: `QRGapFit`, `BlockIncrementalQRGapFit`, `ElementIncrementalQRGapFit`.
  * `transform/`: 2-body (`PairDistanceTransformation`, `CoordinationTransformation`), EAM functions (`FSGenPairFunction`, `CoscutoffPairFunction`, `PolycutoffPairFunction`, `SplinePairTransformation`), 3-body (`Angle3bTransformation`, `CosTransformation`, `CutoffJK3bTransformation`, `Distances3bTransformation`, `MeamTransformation`), aggregators (`TwoBodySum`, `ThreeBodySum`).
  * `sparsification/`: `HistogramUniformSparsifier`, `FromFileSparsifier`.
  * `potentials/`: `ScreenedCoulombPotential` (ZBL), `IsolatedAtomPotential`, `SplinePairPotential`.
  * `splines/`: 1D/3D Cubic B-splines, Hermite, Natural cubic splines.
* **`serialization/`**:
  * Polymorphic HDF5 serialization infrastructure (`SerializationNode`, `SerializationRegistry`, and component serializers).
* **`io/`**:
  * Potential loader (`loadPotential`), QUIP XML converter (`QuipXmlConverter`), tabGAP I/O.
* **`utils/`**:
  * High-level fitting driver (`standardGapFit`), parameter structs (`StandardGapParams`, `StandardTabulationParams`).
* **`pybind/`**:
  * Python C-extension interface (`PyJGAP.cpp`).

### Core Design Idioms: `ValuePtr` & Template Architecture

#### 1. Polymorphic Value Semantics (`ValuePtr<T>`)
In interatomic potential codebases, managing polymorphic types via raw pointers or standard smart pointers (`std::unique_ptr<Base>`) often introduces verbose boilerplate, manual heap management, and a loss of regular value semantics. JGAP solves this using `ValuePtr<T>` (`src/jgap/core/ValuePtr.hpp`):

* **The `Cloneable` Concept**: Constrains `T` to classes providing a `.clone()` method returning `T*` and a virtual destructor:
  ```cpp
  template<typename T>
  concept Cloneable = requires(const T& obj) {
      { obj.clone() } -> std::convertible_to<T*>;
  } && std::has_virtual_destructor_v<T>;
  ```
* **Deep Copying by Default**: Copying a `ValuePtr<T>` automatically clones the underlying object, giving composite potentials (like `GapPotential` and `CompositePotential`) clean value semantics without slicing.
* **Zero-Allocation Moves**: Moving a `ValuePtr<T>` transfers pointer ownership without heap allocations or cloning.
* **Elimination of Allocation Boilerplate**: Derived types can be assigned directly to base pointers without explicit `std::make_unique`:
  ```cpp
  ValuePtr<Potential> pot = ScreenedCoulombPotential(...);
  ```
* **Safe Downcasting**: Provides a type-safe `.as<Derived>()` method wrapping `dynamic_cast<Derived*>(...)`.

#### 2. Template Architecture & Dimensionality
JGAP combines compile-time generic programming for high-performance inner loops with runtime dynamic polymorphism for high-level model assembly:

* **Compile-Time Dimensionality (`size_t Dim`)**:
  Coordinate transformations (`TwoBodyTransformation<Dim>`, `ThreeBodyTransformation<Dim>`), covariance kernels (`Kernel<Dim, DerivDim>`), descriptors (`Descriptor<Dim>`), and aggregators (`NBodyAggregator<Dim>`) are templated on spatial/descriptor dimensionality (`size_t Dim`). This enables fixed-size array allocations (`std::array`), loop unrolling, and SIMD auto-vectorization in tight evaluation loops.
* **Concepts & Inlined Kernel Evaluation**:
  C++20 concepts (such as `CKernelOfDim<Dim>`) constrain kernel implementations. Concrete components like `TwoBodyGapComponent<Dim, TKernel>` accept kernel types as template parameters, allowing the compiler to inline kernel covariance and gradient evaluations directly rather than paying the cost of virtual method dispatches per atomic pair or triplet.
* **Dynamic Boundaries at Component Level**:
  To prevent templates from infecting the entire codebase, top-level interfaces (`GapComponent`, `Potential`, `TabGapComponent`) use dynamic polymorphism. For example, `GapPotential` manages a heterogeneous vector of `ValuePtr<GapComponent>`, allowing models combining different body-orders, dimensions, and kernel formulations to be fitted and evaluated interchangeably.
* **Explicit Template Instantiation**:
  For standard dimensions (`Dim = 1, 2, 3, 4`), serializers and implementations are explicitly instantiated in `.cpp` translation units. This keeps header compilation times manageable while ensuring necessary symbols are exported in the shared library for runtime registries and Python bindings.

### Shared Library Requirement
`jgap_lib` is built strictly as a **shared library (`SHARED`)**. The dynamic serialization registry and polymorphic casting across module boundaries rely on shared-library symbol visibility.

---

## 3. How Serialization Works Overall

`jgap` uses a decoupled, non-intrusive polymorphic serialization system built on top of HDF5 (via HighFive). Potential objects and mathematical components can be saved to and loaded from `.jgap.h5` files without embedding file-format logic inside the core physics classes.

### Architectural Primitives

1. **`SerializationNode`**:
   An abstraction layer over HDF5 groups and datasets. It provides type-safe methods to read and write attributes, strings, scalar datasets, and multidimensional arrays:
   ```cpp
   node.writeAttribute("name", "SquaredExpKernel");
   node.writeDataSet("length_scale", kernel.length_scale);
   ```

2. **`Serialization<TBase>` Interface**:
   Each serializable class hierarchy implements a serializer deriving from `Serialization<TBase>`:
   ```cpp
   class MyPotentialSerialization : public Serialization<Potential> {
   public:
       bool serialize(const ValuePtr<Potential>& obj, SerializationNode& node) const override;
       ValuePtr<Potential> deserialize(const SerializationNode& node) const override;
   };
   ```

3. **`SerializationRegistry<TBase>` & Dynamic Self-Registration**:
   Serializers register themselves dynamically into the registry at library load time using the `REGISTER_SERIALIZATION` macro:
   ```cpp
   REGISTER_SERIALIZATION(MyPotentialSerialization, Potential);
   ```
   When saving a potential:
   ```cpp
   SerializationRegistry<Potential>::serialize(my_pot, "potential.jgap.h5");
   ```
   The registry iterates through all registered serializers for `Potential`. The serializer matching the concrete dynamic type serializes the node. During deserialization, the registry inspects the `"name"` attribute on the HDF5 group to instantiate and populate the corresponding derived object.

---

## 4. How to Add Extensions & New Components

### General Guidelines for New Components & Extensions

When adding new components, cutoffs, kernels, or potential types to `jgap`, follow these general principles:

* **Serialization is Required for objects (Potentials and their parts parametrization) that are intended to be saved for later reuse**:
  If a new component, potential, or parameter object is meant to be saved to or loaded from a `.jgap.h5` potential file, you **must** implement its serialization class inheriting from `Serialization<TBase>` (e.g. `Serialization<GapComponent>` or `Serialization<Potential>`) and register it with `REGISTER_SERIALIZATION(...)` in `src/jgap/serialization/`. Without this, the polymorphic `SerializationRegistry` cannot persist or restore the object.
* **Add and Run Unit Tests (`add/run`)**:
  Whenever adding a new component or extending mathematical primitives, add unit tests in `test/unit/` covering forward evaluation, analytical forces and virials, boundary conditions, and serialization round-trips. Always run the unit test suite (`cmake --workflow --preset dev` or `./build/debug/jgap_tests`) to confirm that all tests pass locally.
* **Validation Suite for Existing Logic**:
  Whenever modifying, refactoring, or optimizing **existing logic** (such as descriptor representations, covariance formulas, solvers, existing component evaluation, or tabulation grids), the automated validation suite **must** be executed (`cmake --workflow --preset validation` or `ctest -L validation`). If numerical outputs or formulas intentionally change, the reference baselines or test tolerances must be adjusted accordingly.

---

### 1. Adding a New GAP Component
1. Inherit from `GapComponent` in `src/jgap/core/potentials/gap/component/`:
   ```cpp
   #include "jgap/core/potentials/gap/component/GapComponent.hpp"

   class CustomGapComponent : public GapComponent {
   public:
       std::optional<AtomicQuantities> covariate(const NeighbourLists& nl) const override;
       Matrix<RowMajor> K_MM() const override;
       size_t nSparsePoints() const override;
       Cutoffs getCutoffs() const override;
       std::set<Species> nonZeroCovarianceFor() const override;
       void tabulate(TabulationData& tables) const override;
       CustomGapComponent* clone() const override;
   };
   ```
2. Implement serialization in `src/jgap/serialization/potentials/gap/component/`:
   - Implement `CustomGapComponentSerialization` deriving from `Serialization<GapComponent>`.
   - Register it with `REGISTER_SERIALIZATION(CustomGapComponentSerialization, GapComponent);`.
3. Add unit tests in `test/unit/potentials/` verifying covariate calculations, energy/force evaluation, and HDF5 serialization round-trips.

### 2. Adding a New Cutoff Function
1. Inherit from `CutoffFunction` in `src/jgap/impl/cutoff/`:
   ```cpp
   #include "jgap/core/cutoff/CutoffFunction.hpp"

   class MyCutoff : public CutoffFunction {
   public:
       MyCutoff(double cutoff, double width);
       double evaluate(double r) const override;
       double derivative(double r) const override;
       double getCutoff() const override;
       MyCutoff* clone() const override;
   };
   ```
2. Create `MyCutoffSerialization.cpp` in `src/jgap/serialization/cutoff/`:
   - Implement `serialize` and `deserialize` storing parameters to HDF5.
   - Register via `REGISTER_SERIALIZATION(MyCutoffSerialization, CutoffFunction);`.
3. Add unit tests in `test/unit/cutoff/` verifying evaluation, derivatives, and serialization round-trips.

### 3. Adding a New Covariance Kernel
1. Inherit from `Kernel<Dim, DerivDim>` in `src/jgap/impl/kernels/`:
   ```cpp
   #include "jgap/core/kernels/Kernel.hpp"
   ```
2. Implement covariance calculation, distance functions, and gradients.
3. Add serialization in `src/jgap/serialization/kernels/` with `REGISTER_SERIALIZATION(..., Kernel<Dim, DerivDim>)`.
4. Add unit tests in `test/unit/kernels/`.

### 4. Adding a New Potential Type
1. Inherit from `Potential` in `src/jgap/impl/potentials/`:
   ```cpp
   #include "jgap/core/potentials/Potential.hpp"

   class CustomPotential : public Potential {
   public:
       AtomicQuantity calculateEnergy(const Atoms& atoms) const override;
       Cutoffs getCutoffs() const override;
       CustomPotential* clone() const override;
       void fillTables(TabulationData& table) const override;
   };
   ```
2. Implement `calculateEnergy` to compute total energy, forces (`AtomicQuantity.forces`), and virials (`AtomicQuantity.virials`).
3. Implement `fillTables` to evaluate the potential over tabulation grids for fast spline export.
4. Implement and register its serializer in `src/jgap/serialization/potentials/` with `REGISTER_SERIALIZATION(CustomPotentialSerialization, Potential)`.
5. Add unit tests in `test/unit/potentials/`.

### 5. Exposing Extensions to Python
In `src/pybind/PyJGAP.cpp`:
1. Bind the new class using `py::class_<T>`.
2. Ensure methods returning C++ objects convert safely to Python-managed instances or NumPy arrays.

---

## 5. Testing & Validation Guide

`jgap` utilizes a two-tier verification strategy.

### Tier 1: Unit Tests (Fast Local Feedback)
* Framework: GoogleTest (`jgap_tests`).
* Location: `test/unit/` (covers math primitives, atomic structures, neighbor lists, cutoffs, splines, and HDF5 round-trips).
* **Execution**:
  ```bash
  # Configure and run via Developer workflow preset:
  cmake --workflow --preset dev

  # Or build and run directly:
  cmake --build build/debug --target jgap_tests
  ./build/debug/jgap_tests

  # Run a targeted test filter:
  ./build/debug/jgap_tests --gtest_filter="TestTabGapPotential.*"
  ```

### Tier 2: Automated Validation Suite (Scientific Accuracy)
* Location: `test/validation/`.
* Verifies numerical accuracy against external references (QUIP XML baselines, tabGAP HDF5 tables, and LAMMPS outputs).
* Includes physical elemental property testing ($a_0, B, C_{ij}, E_{\text{vac}}^f, E_{\text{coh}}$) via ASE + `elastic`.

> [!IMPORTANT]
> **Modifying Existing Logic**:
> Whenever you significantly modify existing logic in `jgap` (e.g. core transformations, covariance kernels, fitting solvers, baseline potentials, or spline tabulation), you should run and/or adjust the validation suite. This ensures physical predictions, reference compatibility with QUIP and LAMMPS, and elemental material property thresholds ($a_0, B, C_{ij}, E_{\text{vac}}^f, E_{\text{coh}}$) remain verified.

* **Execution**:
  1. Extract reference datasets (if not already extracted):
     ```bash
     unzip -q test/resources/validation.zip -d test/resources/
     ```
  2. Run the validation suite workflow preset:
     ```bash
     cmake --workflow --preset validation
     ```
  3. Or run using CTest labels:
     ```bash
     ctest --test-dir build/validation -L validation --output-on-failure
     ```
