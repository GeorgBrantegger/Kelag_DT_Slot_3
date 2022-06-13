import numpy as np 
from math import pi 


def Hagen_Poiseuille(dP,dx,r,vis):
    Q = (pi*r**4)/(8*vis)*dP/dx
    return Q


def get_Q_from_dP(P_above,P_below,dx,constants=[1,1]):
    r       = constants[0]
    vis     = constants[1]


