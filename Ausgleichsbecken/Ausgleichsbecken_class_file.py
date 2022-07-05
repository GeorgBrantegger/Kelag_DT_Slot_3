import numpy as np
# from Ausgleichsbecken_functions import FODE_function, get_h_halfstep, get_p_halfstep

#importing pressure conversion function
import sys
import os
current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)
from functions.pressure_conversion import pressure_conversion

def FODE_function(x, h, alpha, p, rho=1000., g=9.81):
    f = x*abs(x)/h*alpha+g-p/(rho*h)
    return f


class Ausgleichsbecken_class:
# units 
    # make sure that units and print units are the same
    # units are used to label graphs and print units are used to have a bearable format when using pythons print()
    area_unit           = r'$\mathrm{m}^2$'
    area_outflux_unit   = r'$\mathrm{m}^2$'
    flux_unit           = r'$\mathrm{m}^3/\mathrm{s}$'
    level_unit          = 'm'
    pressure_unit       = 'Pa'
    time_unit           = 's'
    volume_unit         = r'$\mathrm{m}^3$'

    area_unit_print           = 'm²'
    area_outflux_unit_print   = 'm²'
    flux_unit_print           = 'm³/s'
    level_unit_print          = 'm'
    pressure_unit_print       = 'Pa'
    time_unit_print           = 's'
    volume_unit_print         = 'm³'

# init
    def __init__(self,area,outflux_area,level_min = 0,level_max = np.inf ,timestep = 1):
        self.area           = area              # base area of the rectangular structure
        self.area_outflux   = outflux_area      # area of the outlet towards the pipeline        
        self.level_min      = level_min         # lowest  allowed water level    
        self.level_max      = level_max         # highest allowed water level
        self.timestep       = timestep          # timestep of the simulation

        # initialize for get_info
        self.level          = "--"
        self.influx         = "--"
        self.outflux        = "--"
        self.volume         = "--"
        

# setter
    def set_volume(self):
        self.volume = self.level*self.area 

    def set_initial_level(self,initial_level):
        self.level = initial_level
        self.set_volume()

    def set_influx(self,influx):
        self.influx = influx

    def set_outflux(self,outflux):
        self.outflux = outflux

# getter
    def get_info(self, full = False):
        new_line = '\n'
        
        if full == True:
            # :<10 pads the self.value to be 10 characters wide
            print_str = (f"The cuboid reservoir has the following attributes: {new_line}" 
                f"----------------------------- {new_line}"
                f"Base area             =       {self.area:<10} {self.area_unit_print} {new_line}"
                f"Outflux area          =       {self.area_outflux:<10} {self.area_outflux_unit_print} {new_line}"
                f"Current level         =       {self.level:<10} {self.level_unit_print}{new_line}"
                f"Critical level low    =       {self.level_min:<10} {self.level_unit_print} {new_line}"
                f"Critical level high   =       {self.level_max:<10} {self.level_unit_print} {new_line}"
                f"Volume in reservoir   =       {self.volume:<10} {self.volume_unit_print} {new_line}"
                f"Current influx        =       {self.influx:<10} {self.flux_unit_print} {new_line}"
                f"Current outflux       =       {self.outflux:<10} {self.flux_unit_print} {new_line}"
                f"Simulation timestep   =       {self.timestep:<10} {self.time_unit_print} {new_line}"
                f"----------------------------- {new_line}")
        else:
            # :<10 pads the self.value to be 10 characters wide
            print_str = (f"The current attributes are: {new_line}" 
                f"----------------------------- {new_line}"
                f"Current level         =       {self.level:<10} {self.level_unit_print}{new_line}"
                f"Volume in reservoir   =       {self.volume:<10} {self.volume_unit_print} {new_line}"
                f"Current influx        =       {self.influx:<10} {self.flux_unit_print} {new_line}"
                f"Current outflux       =       {self.outflux:<10} {self.flux_unit_print} {new_line}"
                f"----------------------------- {new_line}")

        print(print_str)


# methods
    def update_level(self,timestep):
        net_flux = self.influx-self.outflux
        delta_V = net_flux*timestep
        new_level = (self.volume+delta_V)/self.area
        return new_level


    def e_RK_4(self):
        yn = self.outflux/self.area_outflux
        h = self.level
        dt = self.timestep
        p,_ = pressure_conversion(self.pressure,self.pressure_unit,'Pa')
        # update to include p_halfstep
        p_hs,_ = pressure_conversion(self.pressure,self.pressure_unit,'Pa')
        alpha = (self.area_outflux/self.area-1)
        h_hs = self.update_level(dt/2)
        Y1 = yn
        Y2 = yn + dt/2*FODE_function(Y1, h, alpha, self.pressure)
        Y3 = yn + dt/2*FODE_function(Y2, h_hs, alpha, p_hs)
        Y4 = yn + dt*FODE_function(Y3, h_hs, alpha, p_hs)
        ynp1 = yn + dt/6*(FODE_function(Y1, h, alpha, p)+2*FODE_function(Y2, h_hs, alpha, p_hs)+ \
            2*FODE_function(Y3, h_hs, alpha, p_hs)+ FODE_function(Y4, h, alpha, p))

        self.outflux = ynp1*self.area_outflux
