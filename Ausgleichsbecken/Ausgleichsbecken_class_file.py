from logging import exception
import numpy as np

#importing pressure conversion function
import sys
import os
current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)
from functions.pressure_conversion import pressure_conversion

def FODE_function(x,h,A,A_a,p,rho,g):
    # (FODE     ... first order differential equation)
    # based on the outflux formula by Andreas Malcherek
        # https://www.youtube.com/watch?v=8HO2LwqOhqQ
        # adapted for a pressurized pipeline into which the reservoir effuses
            # and flow direction
    # x         ... effusion velocity
    # h         ... level in the reservoir
    # A_a       ... Outflux_Area
    # A         ... Reservoir_Area
    # g         ... gravitational acceleration
    # rho       ... density of the liquid in the reservoir 
    f = x*abs(x)/h*(A_a/A-1)+g-p/(rho*h)
    return f


class Ausgleichsbecken_class:
# units 
    # make sure that units and print units are the same
    # units are used to label graphs and print units are used to have a bearable format when using pythons print()
    area_unit           = r'$\mathrm{m}^2$'
    area_outflux_unit   = r'$\mathrm{m}^2$'
    density_unit        = r'$\mathrm{kg}/\mathrm{m}^3$'
    flux_unit           = r'$\mathrm{m}^3/\mathrm{s}$'
    level_unit          = 'm'
    pressure_unit       = 'Pa'
    time_unit           = 's'
    velocity_unit       = r'$\mathrm{m}/\mathrm{s}$'
    volume_unit         = r'$\mathrm{m}^3$'

    area_unit_print         = 'm²'
    area_outflux_unit_print = 'm²'
    density_unit_print      = 'kg/m³'
    flux_unit_print         = 'm³/s'
    level_unit_print        = 'm'
    pressure_unit_print     = '--' # will be set by .set_pressure() method
    time_unit_print         = 's'
    velocity_unit_print     = 'm/s'
    volume_unit_print       = 'm³'

    g = 9.81    # m/s²  gravitational acceleration


# init
    def __init__(self,area,outflux_area,level_min = 0,level_max = np.inf ,timestep = 1,rho = 1000):
        self.area           = area              # base area of the rectangular structure
        self.area_outflux   = outflux_area      # area of the outlet towards the pipeline        
        self.density        = rho               # density of the liquid in the system
        self.level_min      = level_min         # lowest  allowed water level    
        self.level_max      = level_max         # highest allowed water level
        self.timestep       = timestep          # timestep of the simulation

        # initialize for get_info
        self.influx         = "--"
        self.level          = "--"
        self.outflux        = "--"
        self.volume         = "--"
        

# setter
    def set_initial_level(self,initial_level):
        # sets the level in the reservoir and should only be called during initialization
        if self.level == '--':
            self.level  = initial_level
            self.volume = self.update_volume()
        else:
            raise Exception('Initial level was already set once. Use the .update_level(self,timestep) method to update level based on net flux.')

    def set_influx(self,influx):
        # sets influx to the reservoir in m³/s
            # positive influx means that liquid flows into the reservoir
        self.influx = influx

    def set_outflux(self,outflux):
        # sets outflux to the reservoir in m³/s
            # positive outflux means that liquid flows out of reservoir the reservoir
        self.outflux = outflux

    def set_initial_pressure(self,pressure,display_pressure_unit):
        # sets the static pressure present at the outlet of the reservoir
            # units are used to convert and display the pressure
        self.pressure               = pressure
        self.pressure_unit_print    = display_pressure_unit

    def set_pressure(self,pressure):
        # sets the static pressure present at the outlet of the reservoir
            # units are used to convert and display the pressure
        self.pressure               = pressure

    def set_steady_state(self,ss_influx,ss_level,display_pressure_unit):
        # set the steady state (ss) condition in which the net flux is zero
            # set pressure acting on the outflux area so that the level stays constant
        ss_outflux = ss_influx
        ss_pressure = self.density*self.g*ss_level-(ss_outflux/self.area_outflux)**2*self.density/2

        self.set_influx(ss_influx)
        self.set_initial_level(ss_level)
        self.set_initial_pressure(ss_pressure,display_pressure_unit)
        self.set_outflux(ss_outflux)

# getter
    def get_info(self, full = False):
        new_line = '\n'
        p = pressure_conversion(self.pressure,self.pressure_unit,self.pressure_unit_print)

        
        if full == True:
            # :<10 pads the self.value to be 10 characters wide
            print_str = (f"The cuboid reservoir has the following attributes: {new_line}" 
                f"----------------------------- {new_line}"
                f"Base area             =       {self.area:<10} {self.area_unit_print} {new_line}"
                f"Outflux area          =       {round(self.area_outflux,3):<10} {self.area_outflux_unit_print} {new_line}"
                f"Current level         =       {self.level:<10} {self.level_unit_print}{new_line}"
                f"Critical level low    =       {self.level_min:<10} {self.level_unit_print} {new_line}"
                f"Critical level high   =       {self.level_max:<10} {self.level_unit_print} {new_line}"
                f"Volume in reservoir   =       {self.volume:<10} {self.volume_unit_print} {new_line}"
                f"Current influx        =       {self.influx:<10} {self.flux_unit_print} {new_line}" 
                f"Current outflux       =       {self.outflux:<10} {self.flux_unit_print} {new_line}"
                f"Current outflux vel   =       {round(self.outflux_vel,3):<10} {self.velocity_unit_print} {new_line}"
                f"Current pipe pressure =       {round(p,3):<10} {self.pressure_unit_print} {new_line}"
                f"Simulation timestep   =       {self.timestep:<10} {self.time_unit_print} {new_line}"
                f"Density of liquid     =       {self.density:<10} {self.density_unit_print} {new_line}"
                f"----------------------------- {new_line}")
        else:
            # :<10 pads the self.value to be 10 characters wide
            print_str = (f"The current attributes are: {new_line}" 
                f"----------------------------- {new_line}"
                f"Current level         =       {self.level:<10} {self.level_unit_print}{new_line}"
                f"Volume in reservoir   =       {self.volume:<10} {self.volume_unit_print} {new_line}"
                f"Current influx        =       {self.influx:<10} {self.flux_unit_print} {new_line}"
                f"Current outflux       =       {self.outflux:<10} {self.flux_unit_print} {new_line}"
                f"Current outflux vel   =       {round(self.outflux_vel,3):<10} {self.velocity_unit_print} {new_line}"
                f"Current pipe pressure =       {round(p,3):<10} {self.pressure_unit_print} {new_line}"
                f"----------------------------- {new_line}")

        print(print_str)

    def get_current_level(self):
        return self.level
    
    def get_current_influx(self):
        return self.influx
    
    def get_current_outflux(self):
        return self.outflux
    
    def get_current_volume(self):
        return self.volume
    
    def get_current_pressure(self):
        return self.pressure



# methods
        
    def update_level(self,timestep):
        # update level based on net flux and timestep by calculating the volume change in
            # the timestep and the converting the new volume to a level by assuming a cuboid reservoir
        
        # cannot set new level directly in this method, because it gets called to calcuate during the Runge Kutta
            # to calculate a ficticious level at half the timestep
        net_flux = self.influx-self.outflux
        delta_V = net_flux*timestep
        new_level = (self.volume+delta_V)/self.area
        return new_level

    def update_volume(self):
        # sets volume in reservoir based on self.level
        return self.level*self.area 

    def update_pressure(self):
        p_new = self.density*self.g*self.level-(self.outflux/self.area_outflux)**2*self.density/2
        return p_new

    def timestep_reservoir_evolution(self):
        # update outflux and outflux velocity based on current pipeline pressure and waterlevel in reservoir
        dt      = self.timestep
        rho     = self.density
        g       = self.g
        A       = self.area
        A_a     = self.area_outflux
        yn      = self.outflux/A_a # outflux velocity
        h       = self.level
        h_hs    = self.update_level(dt/2)
        p       = self.pressure
        p_hs    = self.pressure + rho*g*(h_hs-h)
        # explicit 4 step Runge Kutta
        Y1      = yn
        Y2      = yn + dt/2*FODE_function(Y1,h,A,A_a,p,rho,g)
        Y3      = yn + dt/2*FODE_function(Y2,h_hs,A,A_a,p_hs,rho,g)
        Y4      = yn + dt*FODE_function(Y3,h_hs,A,A_a,p_hs,rho,g)
        ynp1    = yn + dt/6*(FODE_function(Y1,h,A,A_a,p,rho,g)+2*FODE_function(Y2,h_hs,A,A_a,p_hs,rho,g)+ \
            2*FODE_function(Y3,h_hs,A,A_a,p_hs,rho,g)+ FODE_function(Y4,h,A,A_a,p,rho,g))

        self.outflux    = ynp1*A_a
        self.level      = self.update_level(dt)
        self.volume     = self.update_volume()
        self.pressure   = self.update_pressure()
        
