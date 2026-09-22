import numpy as np
from ase.calculators.calculator import Calculator, all_changes
import jgap_ase

class JGAPCalculator(Calculator):
    implemented_properties = ['energy', 'forces', 'stress']

    def __init__(self, potential_paths, **kwargs):
        super().__init__(**kwargs)
        if isinstance(potential_paths, str):
            potential_paths = [potential_paths]
        self.calc = jgap_ase.PyJGAPCalculator(potential_paths)

    def calculate(self, atoms=None, properties=['energy', 'forces', 'stress'], system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        
        positions = self.atoms.get_positions()
        chemical_symbols = self.atoms.get_chemical_symbols()
        cell = self.atoms.cell.array.flatten()
        pbc = self.atoms.pbc

        # The C++ calculator returns energy, forces, and virial
        result = self.calc.calculate(chemical_symbols, positions, cell, pbc)
        
        self.results['energy'] = result['energy']
        self.results['forces'] = result['forces']
        
        # ASE expects stress to be -virial / volume (in Voigt order xx, yy, zz, yz, xz, xy)
        volume = self.atoms.get_volume()
        if volume > 0:
            self.results['stress'] = -result['virial'] / volume
        else:
            self.results['stress'] = np.zeros(6)
