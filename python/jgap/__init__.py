"""
jgap: General Atomic Potential framework for machine-learned interatomic potentials.
"""

from typing import List, Optional, Union
import numpy as np

# Import compiled extension
try:
    from ._jgap import (
        Vector3,
        Species,
        Lattice,
        Virials,
        AtomicQuantity,
        Atoms,
        Cutoffs,
        Potential,
        PerConfigTypeSigmas,
        Regularization,
        RegularizationRules,
        PerConfigTypeRegularizationRules,
        SimpleRegularizationRules,
        ScaledRegularizationRules,
        EamPairFunctionType,
        EamMode,
        StandardGapParams,
        standard_gap_fit,
        StandardTabulationParams,
        standard_tabulation,
        load_potential,
        read_atoms,
        write_atoms,
    )
except ImportError as err:
    try:
        from _jgap import (
            Vector3,
            Species,
            Lattice,
            Virials,
            AtomicQuantity,
            Atoms,
            Cutoffs,
            Potential,
            PerConfigTypeSigmas,
            Regularization,
            RegularizationRules,
            PerConfigTypeRegularizationRules,
            SimpleRegularizationRules,
            ScaledRegularizationRules,
            EamPairFunctionType,
            EamMode,
            StandardGapParams,
            standard_gap_fit,
            StandardTabulationParams,
            standard_tabulation,
            load_potential,
            read_atoms,
            write_atoms,
        )
    except ImportError:
        raise ImportError(
            "Failed to load the compiled _jgap C++ extension. "
            "Please ensure that jgap is built and installed with CMake for your active Python version."
        ) from err


def _atoms_to_ase(self):
    """Convert jgap.Atoms to an ase.Atoms object."""
    from ase import Atoms as ASEAtoms

    cell = self.lattice.to_numpy() if self.lattice is not None else None
    pbc = list(self.pbc)
    positions = self.positions
    symbols = self.symbols

    ase_atoms = ASEAtoms(
        symbols=symbols,
        positions=positions,
        cell=cell,
        pbc=pbc,
    )

    if self.energy is not None or self.forces is not None:
        from ase.calculators.singlepoint import SinglePointCalculator

        results = {}
        if self.energy is not None:
            results["energy"] = self.energy
        if self.forces is not None:
            results["forces"] = self.forces
        if self.virials is not None and cell is not None:
            vol = np.abs(np.linalg.det(cell))
            if vol > 0:
                results["stress"] = -self.virials.to_voigt() / vol

        ase_atoms.calc = SinglePointCalculator(ase_atoms, **results)

    return ase_atoms


@classmethod
def _atoms_from_ase(cls, ase_atoms):
    """Create a jgap.Atoms object from an ase.Atoms object."""
    positions = np.ascontiguousarray(ase_atoms.get_positions(), dtype=np.float64)
    symbols = list(ase_atoms.get_chemical_symbols())
    pbc = [bool(p) for p in ase_atoms.pbc]

    cell = ase_atoms.cell.array
    lat = Lattice(cell) if (cell is not None and not np.allclose(cell, 0)) else None

    atoms = cls(positions=positions, symbols=symbols, lattice=lat, pbc=pbc)

    # Check for energy/forces/stress in calculator if present
    if getattr(ase_atoms, "calc", None) is not None:
        res = getattr(ase_atoms.calc, "results", {})
        if "energy" in res:
            atoms.energy = float(res["energy"])
        if "forces" in res:
            atoms.forces = np.ascontiguousarray(res["forces"], dtype=np.float64)

    return atoms


Atoms.to_ase = _atoms_to_ase
Atoms.from_ase = _atoms_from_ase

try:
    from .ase import JGAPCalculator
except ImportError:
    pass

__all__ = [
    "Vector3",
    "Species",
    "Lattice",
    "Virials",
    "AtomicQuantity",
    "Atoms",
    "Cutoffs",
    "Potential",
    "PerConfigTypeSigmas",
    "Regularization",
    "RegularizationRules",
    "PerConfigTypeRegularizationRules",
    "SimpleRegularizationRules",
    "ScaledRegularizationRules",
    "EamPairFunctionType",
    "EamMode",
    "StandardGapParams",
    "standard_gap_fit",
    "StandardTabulationParams",
    "standard_tabulation",
    "load_potential",
    "read_atoms",
    "write_atoms",
    "JGAPCalculator",
]
