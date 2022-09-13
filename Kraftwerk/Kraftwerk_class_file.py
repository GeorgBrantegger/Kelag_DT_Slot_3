import numpy as np
#importing Druckrohrleitung
import sys
import os
current = os.path.dirname(os.path.realpath('Main_Programm.ipynb'))
parent = os.path.dirname(current)
sys.path.append(parent)
from functions.pressure_conversion import pressure_conversion
from Turbinen.Turbinen_class_file import Francis_Turbine

class Kraftwerk_class:
    g = 9.81

    def __init__(self):
        self.turbines = []
        self.n_turbines = 0

# setter
    def set_LAs(self,LA_vec,display_warning=True):
        for i in range(self.n_turbines):
            self.turbines[i].set_LA(LA_vec[i],display_warning)

    def set_pressure(self,pressure):
        for i in range(self.n_turbines):
            self.turbines[i].set_pressure(pressure)
        
    def set_steady_state(self,ss_flux,ss_pressure):
        self.identify_Q_proportion()
        for i in range(self.n_turbines):
            self.turbines[i].set_steady_state(ss_flux*self.Q_prop[i],ss_pressure)

# getter
    def get_current_Q(self):
        Q = 0.
        for i in range(self.n_turbines):
            Q += self.turbines[i].get_current_Q()
        return Q

    def get_current_LAs(self):
        LAs = []
        for i in range(self.n_turbines):
            LAs.append(self.turbines[i].get_current_LA())
        return np.array(LAs)

    def get_current_pressure(self):
        pressures = []
        for i in range(self.n_turbines):
            pressures.append(self.turbines[i].get_current_pressure())
        return np.array(pressures) # consider taking the average, after evaluating how the converge() method affects the result

    def get_n_turbines(self):
        return self.n_turbines

    def get_info(self):
        for turbine in self.turbines:
            turbine.get_info(full=True)

# methods
    def identify_Q_proportion(self):
        Q_n_vec = np.zeros(self.n_turbines)
        for i in range(self.n_turbines):
            Q_n_vec[i] = self.turbines[i].get_Q_n()
        self.Q_prop = Q_n_vec/np.sum(Q_n_vec)

    def add_turbine(self,turbine):
        self.turbines.append(turbine)
        self.n_turbines += 1
    
    def update_LAs(self,LA_soll_vec):
        for i in range(self.n_turbines):
            self.turbines[i].update_LA(LA_soll_vec[i])

    def converge(self,convergence_parameters):
        # small numerical disturbances (~1e-12 m/s) in the velocity can get amplified at the turbine node, because the new velocity of the turbine and the
        # new pressure from the forward characteristic are not perfectly compatible.
        # Therefore, iterate the flux and the pressure so long, until they converge

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
        p_old               = convergence_parameters[9]                     # pressure of previous timestep
        Q_old               = self.get_current_Q()                          
        v_old               = Q_old/area_pipe                               


        while iteration_change > eps:
            
            p_new = p-rho*c*(v_old-v)+rho*c*dt*g*np.sin(alpha)-f_D*rho*c*dt/(2*D)*abs(v)*v
            # print(p_new)
            p_new = p_old+(p_new-p_old)/3
            # print(p_new)
            self.set_pressure(p_new)
            Q_new = self.get_current_Q()
            v_new = Q_new/area_pipe
            # print(Q_old,Q_new)

            iteration_change = abs(Q_old-Q_new)
            Q_old = Q_new.copy()
            v_old = v_new.copy()
            p_old = p_new.copy()
            i = i+1
            if i == 1e6:
                print('did not converge')
                break
        # print(i)





