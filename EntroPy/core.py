# -*- coding: utf-8 -*-
"""
Created on Mon Aug  2 14:34:33 2021

@author: Rafael
"""
# External modules
import numpy as _np
from scipy.integrate import simps as _simps
from scipy.integrate import trapz as _trapz
from scipy.interpolate import interp1d as _interp1d

# Local modules
from .constants import *


def S_from_C(Temp,C,method="t"):
    '''
    Calculates entropy from heat capacity     

    Parameters
    ----------
    Temp: array like
        Temperature data, [K] units
    C : array like
        Heat capacity data, [kb] units
    method : str, optional
        DESCRIPTION. Numerical Method used for the integration.
        Options are the composite (t)rapezoidal or the (s)impson’s rules.
        The default is "t".

    Returns
    -------
    S : array like
        Integrated entropy for each T in [J/mol] units
    '''
    
    S=_np.zeros(len(Temp)) 
    if method == "t":
        for ii in range(0,len(Temp)):
            S[ii] = kb*N_A*_trapz((C[:ii+1])/Temp[:ii+1],x=Temp[:ii+1])
    elif method == "s":
        for ii in range(0,len(Temp)):
            S[ii] = kb*N_A*_simps((C[:ii+1])/Temp[:ii+1],x=Temp[:ii+1])
    return S


def dist_FD (Energy,Temp,mu=0):
    '''
    # Fermi-Dirac distribution f(E,T)

    Parameters
    ----------
    Energy : array like
        Energy data, [kb] units
    Temp : array like
        Temperature data, [K] units
    mu : float, optional
        Total chemical potential. The default is 0.

    Returns
    -------
    f(E,T) : array like
        DESCRIPTION.

    '''
    T=_np.atleast_1d(Temp) #Temperature
    T=_np.where( T > 0, T, _np.nan )
    E=_np.array(Energy)-mu 
    
    #Case T > 0
    
    # Makes use of numerical trick to bypass overflow exp:
    # 1/(1+e^x) == e^(-x)/(1+e^(-x))
    
    condlist = [E > 0,E <= 0]
    
    E_plus=_np.outer(1/T,_np.where(E > 0, E,0))
    E_minus=_np.outer(1/T,_np.where(E <= 0, E,0))
    
    choicelist = [_np.exp(-E_plus)/(_np.exp(-E_plus)+1.0),\
                  1/(_np.exp(E_minus)+1.0)]
    
    result=_np.select(condlist,choicelist)
    
    #Case T = 0
    result[_np.where(_np.isnan(T))]=_np.where(E>0,0,1)
        
    return _np.array(result)


def dist_BE (Energy,Temp,mu=0):
    '''
    Bose-Einstein distribution g(E,T)

    Parameters
    ----------
    Energy : TYPE
        Energy, [kb] units
    Temp : TYPE
        Temperature, [k] units
    mu : TYPE, optional
        Chemical potential in kb units. The default is 0.

    Returns
    -------
    g(E,T) : array like
        DESCRIPTION.

    '''
    T=_np.atleast_1d(Temp)    
    #T=_np.where( _np.abs(T)<1E-12, 1E-12, T )
    #T=_np.where( T >= 1E-12, T, _np.nan )
    
    
    E=_np.array(Energy)-mu
    E=_np.where( _np.abs(E)<1E-12, 1E-12, E)
    E=_np.where( E >= 1E-12, E, _np.nan )

    exp = _np.outer(1/T,E)
    #Case T > 0 and E > 0
    #result= 1/_np.expm1(exp)
    #result= 1/(_np.exp(_np.outer(1/T,E))-1.0)
    
    # 1/(e^x-1) == e^(-x)/(1-e^(-x)) , numerical trick to bypass overflow 
    result= _np.exp(-1*exp)/(-_np.expm1(-1*exp))
    #result= _np.exp(-1*exp)/(1-_np.exp(-1*exp))
    

    return _np.array(result)

def apply_gaussian_smearing(energy,dos, sigma):
    """_summary_

    Args:
        dos (_type_): _description_
        sigma (_type_): _description_

    Returns:
        _type_: _description_
    """

    # Generate the Gaussian filter
    num_points = int(7 * sigma / (energy[1] - energy[0]) + 1)
    x = _np.linspace(-3*sigma, 3*sigma, num_points)
    gaussian_filter = _np.exp(-0.5 * (x / sigma) ** 2) / (sigma * _np.sqrt(2*_np.pi))

    # Perform convolution
    smeared_dos = _np.convolve(dos, gaussian_filter, mode='same')

    # Normalize the smeared DOS
    smeared_dos *= _trapz(dos, energy) / _trapz(smeared_dos, energy)

    return smeared_dos

def S_bosons (Energy,Dos,Temp,units='eV',Ndiv=25) :
    '''
    Entropy for bosonic excitations

    Parameters
    ----------
    Energy : TYPE
        Energy in kb units
    Dos : TYPE
        Energy in kb units
    Temp : TYPE
        DESCRIPTION.
    units : TYPE, optional
        The default is eV.

    Returns
    -------
    TYPE
        DESCRIPTION.

    '''
    T=_np.atleast_1d(Temp)
    
    if units == 'eV' :
        E=Energy/kb_eV
        D=Dos*kb_eV
    elif units == 'Ry' :
        E=Energy/kb_Ry
        D=Dos*kb_Ry
    elif units == 'J' :
        E=Energy/kb
        D=Dos*kb
    elif units == 'THz' :
        E=1E12*h*Energy/kb
        D=Dos*kb/(h*1E12)
    elif units == 'kb' :
        E=Energy
        D=Dos
    
     # Integration of DOS can be tricky, thus the interpolation to improve
    # the S calculation
    Dos_fx=_interp1d(E,D,kind='linear')
    E_new = _np.linspace(E.min(),E.max(),Ndiv*len(Dos))
    Dos_new=Dos_fx(E_new)



    S=Dos_new*((1+dist_BE(E_new,T))*_np.log(1+dist_BE(E_new,T))\
        - dist_BE(E_new,T)*_np.log(dist_BE(E_new,T)))
    return _trapz(S,x=E_new,axis=1)*N_A*kb
    
def S_phonons(Energy,Dos,Temp,units='THz',Ndiv=25):

    return 3*S_bosons(Energy,Dos,Temp,units=units,Ndiv=Ndiv)    

def S_ele (Dos_x,Dos_y,Temp,Ef=0,units='eV',method='t',Ndiv=25):
    '''

    Parameters
    ----------
    Dos_x : TYPE
        DESCRIPTION.
    Dos_y : TYPE
        DESCRIPTION.
    Temp : TYPE
        DESCRIPTION.
    Ef : TYPE, optional
        DESCRIPTION. The default is 0.
    units : TYPE, optional
        DESCRIPTION. The default is 'eV'.
    method : TYPE, optional
        DESCRIPTION. The default is 't'.
    Ndiv : TYPE, optional
        DESCRIPTION. The default is 25.

    Returns
    -------
    S : TYPE
        DESCRIPTION.

    '''
    
    T=_np.array(Temp)
    
    # DOS Conversion to kb units
    E = _np.array(Dos_x) - Ef
    Dos = _np.array(Dos_y)
    
    if units == 'eV' :
        E=E/kb_eV
        Dos=Dos*kb_eV
    elif units == 'Ry' :
        E=E/kb_Ry
        Dos=Dos*kb_Ry
    elif units == 'J' :
        E=E/kb
        Dos=Dos*kb
    elif units == 'kb' :
        pass  
    
    # Integration of DOS can be tricky, thus the interpolation to improve
    # the S calculation
    Dos_fx=_interp1d(E,Dos,kind='linear')
    E_new = _np.linspace(E.min(),E.max(),Ndiv*len(Dos))
    Dos_new=Dos_fx(E_new)
    
    # Calculating mixing term 
    #(1-x)*_np.log(1-x)+x*_np.log(x)
    fd=dist_FD(E_new,T)
    fd=_np.where(_np.logical_or(_np.isclose(fd,0,atol=1e-16),\
                              _np.isclose(fd,1.0,atol=1e-16)), _np.nan,fd)
    mix = _np.where(_np.isfinite(fd),(1-fd)*_np.log(1-fd)+fd*_np.log(fd),0.0)
    
    if method == 't':
        S = -_trapz(mix*Dos_new,x=E_new,axis=1)
    elif method == 's':
        S = -_simps(mix*Dos_new,x=E_new,axis=1)
    return S*N_A*kb
    
def F_ele (Dos_x,Dos_y,Temp,Ef=0,units='eV',method='t',Ndiv=25):
    '''

    Parameters
    ----------
    Dos_x : TYPE
        DESCRIPTION.
    Dos_y : TYPE
        DESCRIPTION.
    Temp : TYPE
        DESCRIPTION.
    Ef : TYPE, optional
        DESCRIPTION. The default is 0.
    units : TYPE, optional
        DESCRIPTION. The default is 'eV'.
    method : TYPE, optional
        DESCRIPTION. The default is 't'.
    Ndiv : TYPE, optional
        DESCRIPTION. The default is 25.

    Returns
    -------
    S : TYPE
        DESCRIPTION.

    '''
    
    T=_np.array(Temp)
    
    # DOS Conversion to kb units
    E = _np.array(Dos_x) - Ef
    Dos = _np.array(Dos_y)
    
    if units == 'eV' :
        E=E/kb_eV
        Dos=Dos*kb_eV
    elif units == 'Ry' :
        E=E/kb_Ry
        Dos=Dos*kb_Ry
    elif units == 'J' :
        E=E/kb
        Dos=Dos*kb
    elif units == 'kb' :
        pass  
    
    # Integration of DOS can be tricky, thus the interpolation to improve
    # the S calculation
    Dos_fx=_interp1d(E,Dos,kind='linear')
    E_new = _np.linspace(E.min(),E.max(),Ndiv*len(Dos))
    Dos_new=Dos_fx(E_new)
    
    # Calculating mixing term 
    #(1-x)*_np.log(1-x)+x*_np.log(x)
    fd=dist_FD(E_new,T)
    
    if method == 't':
        U = _trapz(E_new*fd*Dos_new,x=E_new,axis=1)
        U0 = _trapz(E_new*dist_FD(E_new,0)*Dos_new,x=E_new,axis=1)
        S = S_ele(E_new,Dos_new,T,units='kb',Ndiv=1,method='t')  
    elif method == 's':
        U = _simps(E_new*fd*Dos_new,x=E_new,axis=1)
        U0 = _simps(E_new*dist_FD(E_new,0)*Dos_new,x=E_new,axis=1)
        S = S_ele(E_new,Dos_new,T,units='kb',Ndiv=1,method='s')  
    U-= U0
    return (U-T*S)*N_A*kb

#__all__ = ['apply_gaussian_smearing','dist_BE','dist_FD','S_bosons','S_ele','F_ele','S_from_C']