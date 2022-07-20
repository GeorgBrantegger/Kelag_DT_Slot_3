from matplotlib.pyplot import fill
import numpy as np
from scipy.interpolate import interp2d

#importing pressure conversion function
import sys
import os
current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)
from functions.pressure_conversion import pressure_conversion

class Francis_turbine_class:
    def __init__(self,CSV_name='Durchflusskennlinie.csv'):
        csv = np.genfromtxt(CSV_name,delimiter=',')
        n_rows,_ = np.shape(csv)
        self.raw_csv = np.append(csv,np.zeros([n_rows,1]),axis = 1)
        
    def extract_csv(self,CSV_pressure_unit='bar'):
        ps_vec,_            = pressure_conversion(self.raw_csv[0,1:],CSV_pressure_unit,'Pa')
        self.raw_ps_vec     = np.flip(ps_vec)
        self.raw_LA_vec     = self.raw_csv[1:,0]
        self.raw_Qs_mat     = np.fliplr(self.raw_csv[1:,1:])/1000. # convert from l/s to m³/s

    def get_Q_fun(self):
        Q_fun = interp2d(self.raw_ps_vec,self.raw_LA_vec,self.raw_Qs_mat,bounds_error=False,fill_value=None)
        return Q_fun






