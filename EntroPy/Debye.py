# -*- coding: utf-8 -*-
"""
Created on Mon Aug  2 19:29:22 2021

@author: Rafael
"""
# External modules
import numpy as _np

# Local modules
from .constants import N_A as _N_A
from .constants import kb as _kb

def vdos(v,Debye_T : float): 
    """ Debye Model VDOS 
    
    Returns vibrational DOS according to Debye Model in [states/THz/atom]
    
    Parameters
    ----------
    v : array_like
        Frequencies [THz]
            
    Debye_T : float
              Debye Temperature [K]
    
    Returns
    -------
    dos : array
           Debye VDOS [states/THz/atom]
    """
    
    
    Debye_v = (1E-12)*kb*Debye_T/h # Units [THz]
    
    # If v > 0 and v < Debye_v , VDOS=9*v**2/Debye_v**3
    # Otherwise, VDOS=0
    dos=np.where(np.logical_and(v <= Debye_v, v > 0),\
                  9.0*np.power(v,2)/((Debye_v)**3.0),0.0)
    
    return dos

def S_vib(Temp, Debye_T, N : int = 250):
    '''
    
    Parameters
    ----------
    Temp : float, optional
        Temperature [K]
    Debye_T : float
        Debye Tempereature [K]
    N : int, optional
        Number of sub divisions for integration. The default is 100.

    Returns
    -------
    None.

    '''
    T=np.array(Temp)
    
    A = -3.0*np.log(1-np.exp(-Debye_T/T))
    B = 12.0*np.power(Debye_T/T,-3.0)
    C = []
    for ii in Debye_T/T:
        xx=np.linspace(1E-8,ii,100)
        C.append(np.trapz((xx**3.0)/(np.expm1(xx)),x=xx))
    C = np.array(C) 

    return N_A*kb*(A+B*C)


def F_vib (T,Debye_T):
    # Vibrational Helmoholtz free nergy in Kb/atom units
    
    x = Debye_T/T
    A = (9.0/8.0)*x + 3*np.log(1-np.exp(-x))
    
    Dx = np.linspace(1E-10,x,num=1000)
    Dy = (3.0/np.power(x,3))*np.trapz(np.power(Dx,3)/np.expm1(Dx),x=Dx,axis=0)
    
    return T*(A-Dy)
    