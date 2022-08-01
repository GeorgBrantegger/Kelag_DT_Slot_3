from time import time
import numpy as np
#importing pressure conversion function
import sys
import os
current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)
from functions.pressure_conversion import pressure_conversion

class Francis_Turbine:
 # units 
    # make sure that units and print units are the same
    # units are used to label graphs and print units are used to have a bearable format when using pythons print()
    density_unit        = r'$\mathrm{kg}/\mathrm{m}^3$'
    flux_unit           = r'$\mathrm{m}^3/\mathrm{s}$'
    LA_unit             = '%'
    pressure_unit       = 'Pa'
    time_unit           = 's'
    velocity_unit       = r'$\mathrm{m}/\mathrm{s}$'
    volume_unit         = r'$\mathrm{m}^3$'

    density_unit_print      = 'kg/m³'
    flux_unit_print         = 'm³/s'
    LA_unit_print           = '%'
    pressure_unit_print     = 'mWS' 
    time_unit_print         = 's'
    velocity_unit_print     = 'm/s'
    volume_unit_print       = 'm³'

    g = 9.81    # m/s²  gravitational acceleration

 # init
    def __init__(self, Q_nenn,p_nenn,t_closing=-1.,timestep=-1.):
        self.Q_n    = Q_nenn                                    # nominal flux 
        self.p_n    = p_nenn                                    # nominal pressure
        self.LA_n   = 1. # 100%                                 # nominal Leitapparatöffnung
        h           = pressure_conversion(p_nenn,'Pa','MWs')    # nominal pressure in terms of hydraulic head
        self.A      = Q_nenn/(np.sqrt(2*self.g*h)*0.98)         # Ersatzfläche

        self.dt     = timestep                                  # simulation timestep
        self.t_c = t_closing                                    # closing time
        self.d_LA_max_dt = 1/t_closing                          # maximal change of LA per second

        # initialize for get_info() - parameters will be converted to display -1 if not overwritten
        self.p  = pressure_conversion(-1,self.pressure_unit_print,self.pressure_unit)
        self.Q  = -1.
        self.LA = -0.01 


# setter
    def set_LA(self,LA,display_warning=True):
        # set Leitapparatöffnung
        self.LA = LA
        # warn user, that the .set_LA() method should not be used ot set LA manually
        if display_warning == True:
            print('Consider using the .update_LA() method instead of setting LA manually')

    def set_timestep(self,timestep,display_warning=True):
        # set Leitapparatöffnung
        self.dt = time
        # warn user, that the .set_LA() method should not be used ot set LA manually
        if display_warning == True:
            print('WARNING: You are changing the timestep of the turbine simulation. This has implications on the simulated closing speed!')

    def set_pressure(self,pressure):
        # set pressure in front of the turbine
        self.p = pressure

#getter
    def get_current_Q(self):
        # return the flux through the turbine, based on the current pressure in front
            #  of the turbine and the Leitapparatöffnung
        if self.p < 0:
            self.Q = 0
        else:
            self.Q = self.Q_n*(self.LA/self.LA_n)*np.sqrt(self.p/self.p_n)
        return self.Q

    def get_current_LA(self):
        return self.LA

    def get_info(self, full = False):
        new_line = '\n'
        p   = pressure_conversion(self.p,self.pressure_unit,self.pressure_unit_print)
        p_n = pressure_conversion(self.p_n,self.pressure_unit,self.pressure_unit_print)

        
        if full == True:
            # :<10 pads the self.value to be 10 characters wide
            print_str = (f"Turbine has the following attributes: {new_line}" 
                f"----------------------------- {new_line}"
                f"Type                  =       Francis {new_line}"
                f"Nominal flux          =       {self.Q_n:<10} {self.flux_unit_print} {new_line}"
                f"Nominal pressure      =       {round(p_n,3):<10} {self.pressure_unit_print}{new_line}"
                f"Nominal LA            =       {self.LA_n*100:<10} {self.LA_unit_print} {new_line}"
                f"Closing time          =       {self.t_c:<10} {self.time_unit_print} {new_line}"
                f"Current flux          =       {self.Q:<10} {self.flux_unit_print} {new_line}" 
                f"Current pipe pressure =       {round(p,3):<10} {self.pressure_unit_print} {new_line}"
                f"Current LA            =       {self.LA*100:<10} {self.LA_unit_print} {new_line}"
                f"Simulation timestep   =       {self.dt:<10} {self.time_unit_print} {new_line}"
                f"----------------------------- {new_line}")
        else:
            # :<10 pads the self.value to be 10 characters wide
            print_str = (f"The current attributes are: {new_line}" 
                f"----------------------------- {new_line}"
                f"Current flux          =       {self.Q:<10} {self.flux_unit_print} {new_line}" 
                f"Current pipe pressure =       {round(p,3):<10} {self.pressure_unit_print} {new_line}"
                f"Current LA            =       {self.LA*100:<10} {self.LA_unit_print} {new_line}"
                f"----------------------------- {new_line}")

        print(print_str)

# methods
    def update_LA(self,LA_soll):
        # update the Leitappartöffnung and consider the restrictions of the closing time of the turbine
        LA_diff = self.LA-LA_soll                   # calculate the difference to the target LA
        LA_diff_max = self.d_LA_max_dt*self.dt     # calculate the maximum change in LA based on the given timestep  
        LA_diff = np.sign(LA_diff)*np.min(np.abs([LA_diff,LA_diff_max]))    # calulate the correct change in LA
        
        LA_new = self.LA-LA_diff
        if LA_new < 0.:
            LA_new = 0.
        elif LA_new > 1.:
            LA_new = 1.
        self.set_LA(LA_new,display_warning=False)

    def set_steady_state(self,ss_flux,ss_pressure):
        # calculate and set steady state LA, that allows the flow of ss_flux at ss_pressure through the
            # turbine at the steady state LA
        ss_LA = self.LA_n*ss_flux/self.Q_n*np.sqrt(self.p_n/ss_pressure)
        if ss_LA < 0 or ss_LA > 1:
            raise Exception('LA out of range [0;1]')
        self.set_LA(ss_LA,display_warning=False)
