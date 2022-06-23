import numpy as np 
from math import pi 


def Hagen_Poiseuille(P_above,P_below,dx,constants=[1,1]):
    dP      = P_above-P_below
    r       = constants[0]
    vis     = constants[1]    
    Q = (pi*r**4)/(8*vis)*dP/dx
    return Q


