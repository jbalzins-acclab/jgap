from __future__ import annotations

from typing import List, Union, Sequence, Optional
import numpy as np
from ase.calculators.calculator import Calculator, all_changes

try:
    from ._jgap import Potential, Atoms, load_potential
except ImportError:
    try:
        from _jgap import Potential, Atoms, load_potential
    except ImportError:
        pass


class JGAPCalculator(Calculator):
    """
    ASE Calculator for JGAP potentials.

    Parameters
    ----------
    potential : Union[Potential, str, List[str]]
        A loaded Potential instance, or a file path (e.g. "pot.jgap.h5", "pot.tabgap.h5"),
        or a list of paths (e.g. ["pot.tabgap.h5", "pot.eam.fs"]).
    """

    implemented_properties = ["energy", "forces", "stress"]

    def __init__(self, potential: Union[Potential, str, List[str]], **kwargs):
        super().__init__(**kwargs)

        if isinstance(potential, Potential):
            self.potential = potential
        elif isinstance(potential, str):
            self.potential = load_potential(potential)
        elif isinstance(potential, (list, tuple)):
            self.potential = load_potential(list(potential))
        else:
            raise ValueError(f"Invalid potential input: {potential}")

    def calculate(
        self,
        atoms=None,
        properties=["energy", "forces", "stress"],
        system_changes=all_changes,
    ):
        super().calculate(atoms, properties, system_changes)

        jgap_atoms = Atoms.from_ase(self.atoms)
        result = self.potential.calculate_energy(jgap_atoms)

        self.results["energy"] = result.energy
        self.results["forces"] = result.forces

        # Voigt stress: -virials / volume (in order xx, yy, zz, yz, xz, xy)
        volume = self.atoms.get_volume()
        if volume > 0:
            self.results["stress"] = -result.virials.to_voigt() / volume
        else:
            self.results["stress"] = np.zeros(6)
