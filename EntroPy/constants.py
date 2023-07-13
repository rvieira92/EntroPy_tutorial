from scipy.constants import physical_constants as _physical_constants

from scipy.constants import N_A
N_A = N_A
h = _physical_constants["Planck constant"][0]

Ry2eV=_physical_constants['Rydberg constant times hc in eV'][0]
eV2Ry=1.0/Ry2eV
Ry2J=_physical_constants['Rydberg constant times hc in J'][0]
eV2J=_physical_constants['electron volt-joule relationship'][0] 

from scipy.constants import k as kb
kb_eV = _physical_constants['Boltzmann constant in eV/K'][0]
kb_Ry = kb_eV*eV2Ry

a0 = _physical_constants['atomic unit of length'][0]
au2A = a0*1E10
A2au = 1.0/au2A

#del globals()["_physical_constants"] # No need to import whole dict