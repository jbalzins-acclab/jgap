#!/usr/bin/env python3
"""
Example: fit a standard 2b+3b+EAM GAP on a training set, serialize it, and tabulate it in Python.

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

    args = parser.parse_args()

    total_start = time.time()
    print(f"Fitting on {args.training_file} using ElementIncrementalQRGapFit (RAM limit: {args.ram_limit} GB)")
    training_data = jgap.read_atoms(args.training_file)

    params = jgap.StandardGapParams(
        seed=120,
        screened_coulomb_dataset_file=args.screened_coulomb_dataset_file,
        approx_ram_limit_gb=args.ram_limit,
    )

    # Configure default parameters
    if params.default_eam is not None:
        params.default_eam.eam_pair_function = jgap.EamPairFunctionType.FSGen3
        params.default_eam.eam_mode = jgap.EamMode.Blind

    if params.default_3b is not None:
        params.default_3b.n_sparse = 500
        # Option: use Distances3bTransformation for 3b instead of Angle (default):
        # params.default_3b.transformation_type = jgap.ThreeBodyTransformationType.Distances

    # Default non-species parameters can also be disabled by erasing them (setting to None):
    # params.default_2b = None

    # Species-specific parameters can be provided (multiple components per species allowed):
    # params.species_3b.append(
    #     jgap.StandardGap3bParams(
    #         species=("Fe", "Fe", "Ni"),
    #         transformation_type=jgap.ThreeBodyTransformationType.Distances,
    #         cutoff=3.7,
    #         cutoff_width=0.6,
    #         n_sparse=500,
    #     )
    # )

    rules = jgap.PerConfigTypeRegularizationRules(
        jgap.PerConfigTypeSigmas(0.001, 0.05, 0.1, 0.02)
    )
    sigmas = rules.determine_for_all(training_data)

    potential_file = f"{args.output_prefix}.jgap.h5"
    fit_start = time.time()
    jgap.standard_gap_fit(potential_file, training_data, sigmas, params)
    fit_time = time.time() - fit_start
    print(f"Saved fitted potential to {potential_file}")

    # Tabulate
    tab_start = time.time()
    jgap.standard_tabulation(potential_file, args.output_prefix)
    tab_time = time.time() - tab_start
    print(f"Saved tabulated potential with prefix {args.output_prefix}")

    total_time = time.time() - total_start
    print(f"Fitting execution time:    {fit_time:.2f} s")
    print(f"Tabulation execution time: {tab_time:.2f} s")
    print(f"Total execution time:      {total_time:.2f} s")


if __name__ == "__main__":
    main()
