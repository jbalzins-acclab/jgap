#!/usr/bin/env python3
"""
Example: fit a standard 2b+3b+EAM GAP on a training set, serialize it, and tabulate it in Python.

Usage:
    python standard_fit.py <training.xyz> <output_prefix> [screened_coulomb_dataset_file]
"""

import sys
import time
import jgap


def main():
    if len(sys.argv) < 3 or len(sys.argv) > 4:
        print(f"Usage: {sys.argv[0]} <training.xyz> <output_prefix> [screened_coulomb_dataset_file]")
        sys.exit(1)

    training_file = sys.argv[1]
    output_prefix = sys.argv[2]
    screened_coulomb_dataset_file = sys.argv[3] if len(sys.argv) == 4 else None

    total_start = time.time()
    print(f"Fitting on {training_file}")
    training_data = jgap.read_atoms(training_file)

    params = jgap.StandardGapParams(
        seed=120,
        screened_coulomb_dataset_file=screened_coulomb_dataset_file,
        eam_pair_function=jgap.EamPairFunctionType.FSGen3,
        eam_mode=jgap.EamMode.Blind,
        n_sparse3=500,
        approx_ram_limit_gb=1.0,
        split_sets=[["Fe"], ["Ni"]],
    )

    rules = jgap.PerConfigTypeRegularizationRules(
        jgap.PerConfigTypeSigmas(0.001, 0.05, 0.1, 0.02)
    )
    sigmas = rules.determine_for_all(training_data)

    potential_file = f"{output_prefix}.jgap.h5"
    fit_start = time.time()
    jgap.standard_gap_fit(potential_file, training_data, sigmas, params)
    fit_time = time.time() - fit_start
    print(f"Saved fitted potential to {potential_file}")

    # Tabulate
    tab_start = time.time()
    jgap.standard_tabulation(potential_file, output_prefix)
    tab_time = time.time() - tab_start
    print(f"Saved tabulated potential with prefix {output_prefix}")

    total_time = time.time() - total_start
    print(f"Fitting execution time:    {fit_time:.2f} s")
    print(f"Tabulation execution time: {tab_time:.2f} s")
    print(f"Total execution time:      {total_time:.2f} s")


if __name__ == "__main__":
    main()
