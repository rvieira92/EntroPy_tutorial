# External modules_N_A
import numpy as _np

# Local modules
from .constants import N_A as _N_A
from .constants import eV2J as _ev2J
from .constants import Ry2eV as _Ry2eV
from .constants import kb as _kb

# Expressions as in  https://doi.org/10.1103/physrevb.52.8813

def S_sommerfeld(Temp,D,units='eV'):
    if units == 'eV' :
        D_new = D/_ev2J
    elif units == 'Ry' :
        D_new =  D/Ry2J
    elif units == 'J' :
        D_new = D
    
    T=np.atleast_1d(Temp) #Temperature
    return _N_A*_kb*_kb*np.pi*np.pi*T*D_new/3.0

def E_sommerfeld(Temp,D,units='eV'):
    if units == 'eV' :
        D_new = D/_ev2J
    elif units == 'Ry' :
        D_new =  D/Ry2J
    elif units == 'J' :
        D_new = D

    T=np.atleast_1d(Temp) #Temperature
    return _N_A*_kb*_kb*np.pi*np.pi*T*T*D_new/6.0

def F_sommerfeld(Temp,D,units='eV'):
    if units == 'eV' :
        D_new = D/_ev2J
    elif units == 'Ry' :
        D_new =  D/Ry2J
    elif units == 'J' :
        D_new = D

    T=np.atleast_1d(Temp) #Temperature
    return E_sommerfeld(T,D,units)-T*S_sommerfeld(T,D,units)

def C_sommerfeld(Temp,D,units='eV'):
    if units == 'eV' :
        D_new = D/_ev2J
    elif units == 'Ry' :
        D_new =  D/Ry2J
    elif units == 'J' :
        D_new = D
    
    T=np.atleast_1d(Temp) #Temperature
    #return _N_A*_kb*_kb*np.pi*np.pi*T*D_new/3.0
    return S_sommerfeld(T,D,units)