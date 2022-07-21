import numpy as np
#importing pressure conversion function
import sys
import os
current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)
from functions.pressure_conversion import pressure_conversion

class Francis_Turbine:
    def __init__(self, Q_nenn,p_nenn):
        self.Q_n    = Q_nenn
        self.p_n    = p_nenn
        self.LA_n   = 1. # 100%
        h,_ = pressure_conversion(p_nenn,'Pa','MWs')
        self.A      = Q_nenn/(np.sqrt(2*9.81*h)*0.98)

    def set_LA(self,LA):
        self.LA = LA

    def get_Q(self,p):
        self.Q = self.Q_n*(self.LA/self.LA_n)*np.sqrt(p/self.p_n)
        return self.Q

    def set_closing_time(self,t_closing):
        self.t_c = t_closing
        self.d_LA_max_dt = 1/t_closing

    def change_LA(self,LA_soll,timestep):
        LA_diff = self.LA-LA_soll
        LA_diff_max = self.d_LA_max_dt*timestep
        if abs(LA_diff) > LA_diff_max:
            LA_diff = np.sign(LA_diff)*LA_diff_max
        self.LA = self.LA-LA_diff

    def set_steady_state(self,ss_flux,ss_pressure):
        ss_LA = self.LA_n*ss_flux/self.Q_n*np.sqrt(self.p_n/ss_pressure)
        self.set_LA(ss_LA)
        if ss_LA < 0 or ss_LA > 1:
            print('LA out of range')
