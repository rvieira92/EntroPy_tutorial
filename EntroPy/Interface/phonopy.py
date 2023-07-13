# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 17:29:30 2020

@author: Rafael
"""
import csv
import numpy as _np
def read_thermal_properties(filename):
    
    with open(filename,"r") as input1:
                 
        for line in input1 :
            line = line.strip()
            if "natom:" in line:
                reader=csv.reader([line],delimiter=" ", skipinitialspace=True)
                data = next(reader)
                Natoms = int(data[-1])
            if "zero_point_energy:" in line:
                reader=csv.reader([line],delimiter=" ", skipinitialspace=True)
                data = next(reader)
                E0 = 1000*float(data[-1])     
            if "thermal_properties:" in line:
                break
                
        T=[]
        F=[]
        S=[]
        C=[]

        for line in input1 :   

            if "- temperature:" in line:    
                reader=csv.reader([line],delimiter=" ", skipinitialspace=True)        
                data = next(reader)  
                T.append(float(data[-1]))
            if "free_energy" in line:
                reader=csv.reader([line],delimiter=" ", skipinitialspace=True)
                data = next(reader)
                F.append((1000*float(data[-1]))/Natoms)
            if "entropy" in line:
                reader=csv.reader([line],delimiter=" ", skipinitialspace=True)
                data = next(reader)
                S.append(float(data[-1])/Natoms)
            if "heat_capacity" in line:
                reader=csv.reader([line],delimiter=" ", skipinitialspace=True)
                data = next(reader)
                C.append(float(data[-1])/Natoms)

    return _np.c_[T,F,S,C]
    
    
    