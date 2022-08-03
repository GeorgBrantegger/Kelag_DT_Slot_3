import numpy as np

#importing pressure conversion function
import sys
import os

from pyparsing import alphanums
current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)
from functions.pressure_conversion import pressure_conversion

class Francis_Turbine:
 # units 
    # make sure that units and print units are the same
    # units are used to label graphs and disp units are used to have a bearable format when using pythons print()
    density_unit        = r'$\mathrm{kg}/\mathrm{m}^3$'
    flux_unit           = r'$\mathrm{m}^3/\mathrm{s}$'
    LA_unit             = '%'
    pressure_unit       = 'Pa'
    time_unit           = 's'
    velocity_unit       = r'$\mathrm{m}/\mathrm{s}$'
    volume_unit         = r'$\mathrm{m}^3$'

    density_unit_disp   = 'kg/m³'
    flux_unit_disp      = 'm³/s'
    LA_unit_disp        = '%'
    time_unit_disp      = 's'
    velocity_unit_disp  = 'm/s'
    volume_unit_disp    = 'm³'

    g = 9.81    # m/s²  gravitational acceleration

 # init
    def __init__(self, Q_nenn,p_nenn,t_closing,timestep,pressure_unit_disp):
        self.Q_n    = Q_nenn                                    # nominal flux 
        self.p_n    = p_nenn                                    # nominal pressure
        self.LA_n   = 1. # 100%                                 # nominal Leitapparatöffnung
        self.dt     = timestep                                  # simulation timestep
        self.t_c    = t_closing                                 # closing time
        self.d_LA_max_dt = 1/t_closing                          # maximal change of LA per second

        self.pressure_unit_disp = pressure_unit_disp

        # initialize for get_info() - parameters will be converted to display -1 if not overwritten
        self.p  = pressure_conversion(-1,self.pressure_unit_disp,self.pressure_unit)
        self.Q  = -1.
        self.LA = -0.01 


# setter
    def set_LA(self,LA,display_warning=True):
        # set Leitapparatöffnung
        self.LA = LA
        # warn user, that the .set_LA() method should not be used ot set LA manually
        if display_warning == True:
            print('You are setting the guide vane opening of the turbine manually. \n \
                This is not an intended use of this method. \n \
                Refer to the .update_LA() method instead.')

    def set_pressure(self,pressure):
        # set pressure in front of the turbine
        self.p = pressure

    def set_steady_state(self,ss_flux,ss_pressure):
        # calculate and set steady state LA, that allows the flow of ss_flux at ss_pressure through the
            # turbine at the steady state LA
        ss_LA = self.LA_n*ss_flux/self.Q_n*np.sqrt(self.p_n/ss_pressure)
        if ss_LA < 0 or ss_LA > 1:
            raise Exception('LA out of range [0;1]')
        self.set_LA(ss_LA,display_warning=False)

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

    def get_current_pressure(self):
        return pressure_conversion(self.p,self.pressure_unit,self.pressure_unit_disp)

    def get_info(self, full = False):
        new_line = '\n'
        p   = pressure_conversion(self.p,self.pressure_unit,self.pressure_unit_disp)
        p_n = pressure_conversion(self.p_n,self.pressure_unit,self.pressure_unit_disp)

        
        if full == True:
            # :<10 pads the self.value to be 10 characters wide
            print_str = (f"Turbine has the following attributes: {new_line}" 
                f"----------------------------- {new_line}"
                f"Type                  =       Francis {new_line}"
                f"Nominal flux          =       {self.Q_n:<10} {self.flux_unit_disp} {new_line}"
                f"Nominal pressure      =       {round(p_n,3):<10} {self.pressure_unit_disp}{new_line}"
                f"Nominal LA            =       {self.LA_n*100:<10} {self.LA_unit_disp} {new_line}"
                f"Closing time          =       {self.t_c:<10} {self.time_unit_disp} {new_line}"
                f"Current flux          =       {self.Q:<10} {self.flux_unit_disp} {new_line}" 
                f"Current pipe pressure =       {round(p,3):<10} {self.pressure_unit_disp} {new_line}"
                f"Current LA            =       {self.LA*100:<10} {self.LA_unit_disp} {new_line}"
                f"Simulation timestep   =       {self.dt:<10} {self.time_unit_disp} {new_line}"
                f"----------------------------- {new_line}")
        else:
            # :<10 pads the self.value to be 10 characters wide
            print_str = (f"The current attributes are: {new_line}" 
                f"----------------------------- {new_line}"
                f"Current flux          =       {self.Q:<10} {self.flux_unit_disp} {new_line}" 
                f"Current pipe pressure =       {round(p,3):<10} {self.pressure_unit_disp} {new_line}"
                f"Current LA            =       {self.LA*100:<10} {self.LA_unit_disp} {new_line}"
                f"----------------------------- {new_line}")

        print(print_str)

# update methods
    def update_LA(self,LA_soll):
        # update the Leitappartöffnung and consider the restrictions of the closing time of the turbine
        LA_diff = self.LA-LA_soll                                           # calculate the difference to the target LA
        LA_diff_max = self.d_LA_max_dt*self.dt                              # calculate the maximum possible change in LA based on the given timestep  
        LA_diff = np.sign(LA_diff)*np.min(np.abs([LA_diff,LA_diff_max]))    # calulate the correct change in LA
        
        # make sure that the LA is not out of the range [0;1]
        LA_new = self.LA-LA_diff
        if LA_new < 0.:
            LA_new = 0.
        elif LA_new > 1.:
            LA_new = 1.
        self.set_LA(LA_new,display_warning=False)

# methods
    def converge(self,convergence_parameters):
        # small numerical disturbances (~1e-12 m/s) in the velocity can get amplified at the turbine node, because the new velocity of the turbine and the
        # new pressure from the forward characteristic are not compatible.
        eps                 = 1e-12                                         # convergence criterion: iteration change < eps
        iteration_change    = 1.                                            # change in Q from one iteration to the next
        i                   = 0                                             # safety variable. break loop if it exceeds 1e6 iterations
        g                   = self.g                                        # gravitational acceleration
        p                   = convergence_parameters[0]                     # pressure at second to last node (see method of characterisctics - boundary condidtions)
        v                   = convergence_parameters[1]                     # velocity at second to last node (see method of characterisctics - boundary condidtions)
        D                   = convergence_parameters[2]                     # diameter of the pipeline
        area_pipe           = convergence_parameters[3]                     # area of the pipeline
        alpha               = convergence_parameters[4]                     # elevation angle of the pipeline
        f_D                 = convergence_parameters[5]                     # Darcy friction coefficient
        c                   = convergence_parameters[6]                     # pressure wave propagtation velocity
        rho                 = convergence_parameters[7]                     # density of the liquid
        dt                  = convergence_parameters[8]                     # timestep of the characteristic method
        
        p_old               = self.get_current_pressure()                   
        Q_old               = self.get_current_Q()                          
        v_old               = Q_old/area_pipe                               


        while iteration_change > eps:
            self.set_pressure(p_old)
            Q_new = self.get_current_Q()
            v_new = Q_new/area_pipe
            p_new = p-rho*c*(v_old-v)+rho*c*dt*g*np.sin(alpha)-f_D*rho*c*dt/(2*D)*abs(v)*v

            iteration_change = abs(Q_old-Q_new)
            Q_old = Q_new.copy()
            p_old = p_new.copy()
            v_old = v_new.copy()
            i = i+1
            if i == 1e6:
                print('did not converge')
                break
        # print(i)
        self.Q = Q_new