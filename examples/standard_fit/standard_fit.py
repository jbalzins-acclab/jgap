#!/usr/bin/env python3
"""
Example: fit a standard 2b+3b+EAM GAP on a training set, serialize it, and tabulate it in Python.

Demonstrates full configuration flexibility:
  1. Custom XYZ property names (energy, force, virial/stress) via jgap.read_atoms.
  2. Hyperparameters for 2-body, EAM, and 3-body (cutoffs, n_sparse, energy_scale, length_scale).
  3. Transformation type selection (Angle vs Distances for 3-body).
  4. Disabling specific components or adding species-specific overrides.
  5. Multi-level regularization with isotropic and anisotropic virials per config_type.
  6. High-resolution tabulation grid tuning via StandardTabulationParams.
  7. Verification: loading potential and predicting energy/forces/virials.

Usage:
    python standard_fit.py <training.xyz> <output_prefix> [screened_coulomb_dataset_file] [--ram-limit <gb>]
"""

import argparse
import sys
import time
import jgap


def main():
    parser = argparse.ArgumentParser(
        description="Fit a standard 2b+3b+EAM GAP potential with ElementIncrementalQRGapFit and tabulate."
    )
    parser.add_argument("training_file", help="Path to training .xyz file")
    parser.add_argument("output_prefix", help="Output prefix for fitted potential")
    parser.add_argument(
        "screened_coulomb_dataset_file",
        nargs="?",
        default=None,
        help="Optional path to screened coulomb dataset file",
    )
    parser.add_argument(
        "--ram-limit",
        type=float,
        default=2.0,
        help="RAM limit in GB for ElementIncrementalQRGapFit out-of-core execution (default: 2.0)",
    )
    parser.add_argument(
        "--virial-prop",
        default="virial",
        help="Property name for virials in extended XYZ (default: 'virial')",
    )
    parser.add_argument(
        "--energy-prop",
        default="energy",
        help="Property name for energy in extended XYZ (default: 'energy')",
    )
    parser.add_argument(
        "--force-prop",
        default="force",
        help="Property name for atomic forces in extended XYZ (default: 'force')",
    )

    args = parser.parse_args()

    total_start = time.time()
    print(f"Reading training data from {args.training_file}...")

    # -------------------------------------------------------------------------
    # 1. Dataset Loading with Custom Property Names
    # -------------------------------------------------------------------------
    # jgap.read_atoms supports overriding any extended XYZ property names:
    # positions, species, forces/force, virials/virial, energy, lattice, pbc, config_type
    training_data = jgap.read_atoms(
        args.training_file,
        virial=args.virial_prop,
        energy=args.energy_prop,
        force=args.force_prop,
    )
    print(f"Loaded {len(training_data)} frames (RAM limit: {args.ram_limit} GB)")

    # -------------------------------------------------------------------------
    # 2. Hyperparameters & Component Configuration (StandardGapParams)
    # -------------------------------------------------------------------------
    params = jgap.StandardGapParams(
        seed=120,
        screened_coulomb_dataset_file=args.screened_coulomb_dataset_file,
        approx_ram_limit_gb=args.ram_limit,
    )

    # --- 2-Body Component (Distance) ---
    params.default_2b.cutoff = 5.2
    params.default_2b.cutoff_width = 1.0
    params.default_2b.n_sparse = 20
    params.default_2b.energy_scale = 10.0      # delta / signal variance (eV)
    params.default_2b.length_scale = 1.0      # radial distance scale (Angstrom)

    # --- EAM Component (Electron Density) ---
    params.default_eam.cutoff = 5.2
    params.default_eam.n_sparse = 20
    params.default_eam.min_density = 0.05
    params.default_eam.eam_pair_function = jgap.EamPairFunctionType.FSGen3
    params.default_eam.eam_mode = jgap.EamMode.Blind  # Blind, FSsym, FSgen, EAM
    params.default_eam.energy_scale = 1.0     # delta (eV)
    params.default_eam.length_scale = 1.0     # density correlation scale

    # --- 3-Body Component (Angle or Distances) ---
    params.default_3b.cutoff = 4.0
    params.default_3b.cutoff_width = 0.6
    params.default_3b.n_sparse = 500
    params.default_3b.energy_scale = 1.0      # delta (eV)
    params.default_3b.length_scale = 1.0      # isotropic length scale (or length_scales=[1.0, 1.0, 1.0])
    # Transformation type: Angle (theta_ijk) or Distances (rij, rik, rjk)
    params.default_3b.transformation_type = jgap.ThreeBodyTransformationType.Angle

    # --- Component Disabling Example ---
    # To fit without 2-body or without EAM, set the default to None:
    # params.default_2b = None

    # --- Species-Specific Custom Overrides (2-Body, EAM, 3-Body) ---
    # When a custom species parameter is added to species_2b, species_eam, or species_3b,
    # that exact species combination is AUTOMATICALLY EXCLUDED from default expansion.
    # The default parameters (e.g. default_2b / default_eam / default_3b) will then only
    # generate components for any remaining species combinations found in the training data.
    # (If default_* is set to None, only the explicitly added species components are fitted.)
    #
    # 1. Custom 2-Body pair override (e.g., Fe-Ni pair with dedicated cutoff and scales):
    # params.species_2b.append(
    #     jgap.StandardGap2bParams(
    #         species=("Fe", "Ni"),
    #         cutoff=4.8,
    #         cutoff_width=0.8,
    #         n_sparse=30,
    #         energy_scale=8.0,
    #         length_scale=1.1,
    #     )
    # )
    #
    # 2. Custom EAM central species override (e.g., dedicated parameters for Ni):
    # params.species_eam.append(
    #     jgap.StandardGapEamParams(
    #         species="Ni",
    #         cutoff=5.0,
    #         n_sparse=25,
    #         min_density=0.08,
    #         eam_pair_function=jgap.EamPairFunctionType.FSGen3,
    #         eam_mode=jgap.EamMode.Blind,
    #         energy_scale=1.2,
    #         length_scale=0.9,
    #     )
    # )
    #
    # 3. Custom 3-Body triplet override (e.g., Fe-Fe-Ni with distance-based transformation):
    # params.species_3b.append(
    #     jgap.StandardGap3bParams(
    #         species=("Fe", "Fe", "Ni"),
    #         transformation_type=jgap.ThreeBodyTransformationType.Distances,
    #         cutoff=3.8,
    #         cutoff_width=0.6,
    #         n_sparse=600,
    #         energy_scale=1.5,
    #         length_scales=[1.2, 1.2, 1.2], # or length_scale=1.2
    #     )
    # )

    # -------------------------------------------------------------------------
    # 3. Regularization & Noise Weights (PerConfigTypeRegularizationRules)
    # -------------------------------------------------------------------------
    # Sigmas are in physical units:
    #   energy:        eV / atom
    #   force:         eV / Angstrom
    #   virials_iso:   isotropic stress / virial (eV / atom)
    #   virials_aniso: deviatoric / anisotropic virial (eV / atom)
    default_sigmas = jgap.PerConfigTypeSigmas(
        energy=0.001,
        force=0.05,
        virials_iso=0.1,
        virials_aniso=0.02,
    )

    # Config_type-specific weights:
    rules = jgap.PerConfigTypeRegularizationRules(
        default_sigmas=default_sigmas,
        exact_config_type_sigmas={
            "isolated_atom": jgap.PerConfigTypeSigmas(energy=0.0001, force=0.01, virials=0.1),
            "IsolatedAtom": jgap.PerConfigTypeSigmas(energy=0.0001, force=0.01, virials=0.1),
        },
        config_type_contains_sigmas={
            "liquid": jgap.PerConfigTypeSigmas(energy=0.005, force=0.1, virials_iso=0.2, virials_aniso=0.05),
            "Liquid": jgap.PerConfigTypeSigmas(energy=0.005, force=0.1, virials_iso=0.2, virials_aniso=0.05),
        },
    )

    # Or via string (config:energy:force:virial:dummy, where dummy added to allow compatibility with the old QUIP format):
    # config_str = "isolated_atom:0.0001:0.01:0.1:0.0:liquid:0.005:0.1:0.2:0.0"
    # rules = jgap.PerConfigTypeRegularizationRules(default_sigmas, config_str)
    sigmas = rules.determine_for_all(training_data)

    # -------------------------------------------------------------------------
    # 4. Fitting Linear GAP (standard_gap_fit)
    # -------------------------------------------------------------------------
    potential_file = f"{args.output_prefix}.jgap.h5"
    print(f"Fitting potential -> {potential_file}...")
    fit_start = time.time()
    jgap.standard_gap_fit(potential_file, training_data, sigmas, params)
    fit_time = time.time() - fit_start
    print(f"Saved fitted potential to {potential_file}")

    # -------------------------------------------------------------------------
    # 5. Tabulation into 3D B-Splines & EAM (standard_tabulation)
    # -------------------------------------------------------------------------
    tab_params = jgap.StandardTabulationParams(
        r_min_3b=0.5,                 # Minimum 3b radius (Angstrom)
        max_eam_density=10.0,         # Maximum EAM density bound
        n_grid_2b=5000,               # 2-body spline grid points
        n_grid_3b=[80, 80, 80],       # 3-body 3D cubic B-spline grid resolution
    )

    print(f"Tabulating potential -> {args.output_prefix}.tabgap.h5 & .eam.fs...")
    tab_start = time.time()
    jgap.standard_tabulation(potential_file, args.output_prefix, tab_params)
    tab_time = time.time() - tab_start
    print(f"Saved tabulated potential with prefix {args.output_prefix}")

    # -------------------------------------------------------------------------
    # 6. Verification & Evaluation
    # -------------------------------------------------------------------------
    print("Verifying fitted & tabulated potentials on test frame...")
    pot_fit = jgap.load_potential(potential_file)
    pot_tab = jgap.load_potential(f"{args.output_prefix}.tabgap.h5")

    test_frame = training_data[0]
    res_fit = pot_fit.calculate_energy(test_frame)
    res_tab = pot_tab.calculate_energy(test_frame)

    print(f"  Fitted GAP energy:    {res_fit.energy:.6f} eV")
    print(f"  Tabulated GAP energy: {res_tab.energy:.6f} eV")
    print(f"  Difference:           {abs(res_fit.energy - res_tab.energy):.6e} eV")

    total_time = time.time() - total_start
    print("\nTiming Summary:")
    print(f"  Fitting time:    {fit_time:.2f} s")
    print(f"  Tabulation time: {tab_time:.2f} s")
    print(f"  Total runtime:   {total_time:.2f} s")


if __name__ == "__main__":
    main()
