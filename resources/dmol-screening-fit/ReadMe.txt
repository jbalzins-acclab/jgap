This folder contains ScreenedCoulombPotential coefficients, based on energy datasets /dmol and /mp2 from:
    https://zenodo.org/records/17302337
("Repulsive interatomic potentials calculated at three levels of theory" by K. Nordlund, S. Lehtola and G. Hobler)

Fit was performed via the script scripts/screened_coulomb/fit_screened_coulomb.py,
which is a modification of a script written by Jesper Byggmästar
that used to be a part of a workflow for creating quip-compatible screened Coulomb tabulated pair-potentials.

Files are written in format:
element_symbol1 element_symbol2 c1 c2 c3 c4 c5 c6 - per line
    ^(corresponding to c1->a1 c2->b1 c3->a2 c4->b2 c5->a3 c6->b3 in the notation of abovementioned paper;
       b_i's in the units of 1/Å)