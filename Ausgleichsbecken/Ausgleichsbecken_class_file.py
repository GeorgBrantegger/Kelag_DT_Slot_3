import os
import sys
from logging import exception

import numpy as np

#importing pressure conversion function
current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)
from functions.pressure_conversion import pressure_conversion


def FODE_function(x_out,h,A,A_a,p,rho,g):
    # (FODE     ... first order differential equation)
    # based on the outflux formula by Andreas Malcherek
        # https://www.youtube.com/watch?v=8HO2LwqOhqQ
        # adapted for a pressurized pipeline into which the reservoir effuses
        #   and flow direction
    # x_out     ... effusion velocity
    # h         ... level in the reservoir
    # A_a       ... Area_outflux
    # A         ... Area_reservoir_base
    # g         ... gravitational acceleration
    # rho       ... density of the liquid in the reservoir 
    f = x_out*abs(x_out)/h*(A_a/A-1.)+g-p/(rho*h)
    return f


class Ausgleichsbecken_class:
# units 
    # make sure that units and display units are the same
    # units are used to label graphs and disp units are used to have a bearable format when using pythons print()
    area_unit               = r'$\mathrm{m}^2$'
    area_outflux_unit       = r'$\mathrm{m}^2$'
    density_unit            = r'$\mathrm{kg}/\mathrm{m}^3$'
    flux_unit               = r'$\mathrm{m}^3/\mathrm{s}$'
    level_unit              = 'm'
    pressure_unit           = 'Pa' # DONT CHANGE needed for pressure conversion
    time_unit               = 's'
    velocity_unit           = r'$\mathrm{m}/\mathrm{s}$'
    volume_unit             = r'$\mathrm{m}^3$'

    area_unit_disp          = 'm²'
    area_outflux_unit_disp  = 'm²'
    density_unit_disp       = 'kg/m³'
    flux_unit_disp          = 'm³/s'
    level_unit_disp         = 'm'
    # pressure_unit_disp will be set within the __init__() method
    time_unit_disp          = 's'
    velocity_unit_disp      = 'm/s'
    volume_unit_disp        = 'm³'

    g = 9.81    # m/s²  gravitational acceleration


# init
    def __init__(self,area,area_outflux,timestep,pressure_unit_disp,level_min=0,level_max=np.inf,rho = 1000.):
        """
        Creates a reservoir with given attributes in this order: \n
            Base Area                       [m²]        \n
            Outflux Area                    [m²]        \n
            Simulation timestep             [s]         \n
            Pressure unit for displaying    [string]    \n
            Minimal level                   [m]         \n
            Maximal level                   [m]         \n
            Density of the liquid           [kg/m³]     \n
        """
        #set initial attributes
        self.area                   = area                  # base area of the cuboid reservoir
        self.area_out               = area_outflux          # area of the outlet towards the pipeline
        self.density                = rho                   # density of the liquid in the system
        self.level_min              = level_min             # lowest  allowed water level - warning yet to be implemented
        self.level_max              = level_max             # highest allowed water level - warning yet to be implemented
        self.pressure_unit_disp     = pressure_unit_disp    # pressure unit for displaying
        self.timestep               = timestep              # timestep in the time evolution method

        # initialize for get_info() (if get_info() gets called before set_steady_state() is executed)
        self.influx     = -np.inf
        self.outflux    = -np.inf
        self.level      = -np.inf
        self.pressure   = -np.inf
        self.volume     = -np.inf
        

# setter - set attributes
    def set_initial_level(self,initial_level):
        # sets the initial level in the reservoir and should only be called during initialization
        if self.level == -np.inf:
            self.level = initial_level
            self.update_volume(set_flag=True)
        else:
            raise Exception('Initial level was already set once. Use the .update_level(self,timestep,set_flag=True) method to update level based on net flux.')

    def set_initial_pressure(self,initial_pressure):
        # sets the initial static pressure present at the outlet of the reservoir and should only be called during initialization
        if self.pressure == -np.inf:
            self.pressure = initial_pressure
        else:
            raise Exception('Initial pressure was already set once. Use the .update_pressure(self) method to update pressure based current level.')

    def set_influx(self,influx):
        # sets influx to the reservoir in m³/s
            # positive influx means that liquid flows into the reservoir
        self.influx = influx

    def set_outflux(self,outflux,display_warning=True):
        # sets outflux to the reservoir in m³/s
            # positive outflux means that liquid flows out of reservoir the reservoir
        if display_warning == True:
            print('You are setting the outflux from the reservoir manually. \n \
                This is not an intended use of this method. \n \
                Refer to the timestep_reservoir_evolution() or set_steady_state() method instead.')
        self.outflux = outflux

    def set_level(self,level,display_warning=True):
        # sets level in the reservoir in m
        if display_warning == True:
            print('You are setting the level of the reservoir manually. \n \
                This is not an intended use of this method. \n \
                Refer to the update_level() or set_steady_state() method instead.')
        self.level = level 

    def set_pressure(self,pressure,display_warning=True):
        # sets pressure in the pipeline just below the reservoir in Pa
        if display_warning == True:
            print('You are setting the pressure below the reservoir manually. \n \
                This is not an intended use of this method. \n \
                Refer to the update_pressure() or set_steady_state() method instead.')
        self.pressure = pressure 

    def set_volume(self,volume,display_warning=True):
        # sets volume in reservoir
        if display_warning == True:
            print('You are setting the volume in the reservoir manually. \n \
                This is not an intended use of this method. \n \
                Refer to the .update_volume() or set_initial_level() or set_steady_state() method instead.')
        self.volume = volume

    def set_steady_state(self,ss_influx,ss_level):
        # set the reservoir to steady state (ss) condition in which the net flux is zero
            # set pressure acting on the outflux area so that the level stays constant
        ss_outflux      = ss_influx
        ss_influx_vel   = abs(ss_influx/self.area)
        ss_outflux_vel  = abs(ss_outflux/self.area_out)
        # see confluence doc for explaination on how to arrive at the ss pressure formula
        ss_pressure = self.density*self.g*ss_level+self.density*ss_outflux_vel*(ss_influx_vel-ss_outflux_vel)

        # use setter methods to set the attributes to their steady state values
        self.set_influx(ss_influx)
        self.set_initial_level(ss_level)
        self.set_initial_pressure(ss_pressure)
        self.set_outflux(ss_outflux,display_warning=False)

# getter - return attributes
    def get_info(self, full = False):
        # prints out the info on the current state of the reservoir
        new_line = '\n'
        if self.pressure != np.inf:
            p = pressure_conversion(self.pressure,self.pressure_unit,self.pressure_unit_disp)
        if self.outflux != np.inf:
            outflux_vel = self.outflux/self.area_out
        
        if full == True:
            # :<10 pads the self.value to be 10 characters wide
            print_str = (f"The cuboid reservoir has the following attributes: {new_line}" 
                f"----------------------------- {new_line}"
                f"Base area             =       {self.area:<10} {self.area_unit_disp} {new_line}"
                f"Outflux area          =       {round(self.area_out,3):<10} {self.area_outflux_unit_disp} {new_line}"
                f"Current level         =       {self.level:<10} {self.level_unit_disp}{new_line}"
                f"Critical level low    =       {self.level_min:<10} {self.level_unit_disp} {new_line}"
                f"Critical level high   =       {self.level_max:<10} {self.level_unit_disp} {new_line}"
                f"Volume in reservoir   =       {self.volume:<10} {self.volume_unit_disp} {new_line}"
                f"Current influx        =       {round(self.influx,3):<10} {self.flux_unit_disp} {new_line}" 
                f"Current outflux       =       {round(self.outflux,3):<10} {self.flux_unit_disp} {new_line}"
                f"Current outflux vel   =       {round(outflux_vel,3):<10} {self.velocity_unit_disp} {new_line}"
                f"Current pipe pressure =       {round(p,3):<10} {self.pressure_unit_disp} {new_line}"
                f"Simulation timestep   =       {self.timestep:<10} {self.time_unit_disp} {new_line}"
                f"Density of liquid     =       {self.density:<10} {self.density_unit_disp} {new_line}"
                f"----------------------------- {new_line}")
        else:
            # :<10 pads the self.value to be 10 characters wide
            print_str = (f"The current attributes are: {new_line}" 
                f"----------------------------- {new_line}"
                f"Current level         =       {self.level:<10} {self.level_unit_disp}{new_line}"
                f"Current volume        =       {self.volume:<10} {self.volume_unit_disp} {new_line}"
                f"Current influx        =       {round(self.influx,3):<10} {self.flux_unit_disp} {new_line}"
                f"Current outflux       =       {round(self.outflux,3):<10} {self.flux_unit_disp} {new_line}"
                f"Current outflux vel   =       {round(outflux_vel,3):<10} {self.velocity_unit_disp} {new_line}"
                f"Current pipe pressure =       {round(p,3):<10} {self.pressure_unit_disp} {new_line}"
                f"----------------------------- {new_line}")

        print(print_str)
    
    def get_current_influx(self):
        return self.influx
    
    def get_current_outflux(self):
        return self.outflux

    def get_current_level(self):
        return self.level
    
    def get_current_pressure(self):
        return self.pressure

    def get_current_volume(self):
        return self.volume

# update methods - update attributes based on some parameter     
    def update_level(self,timestep,set_flag=False):
        # update level based on net flux and timestep by calculating the volume change in
            # the timestep and the converting the new volume to a level by assuming a cuboid reservoir
        net_flux = self.influx-self.outflux
        delta_level = net_flux*timestep/self.area
        level_new = (self.level+delta_level)
        # set flag is necessary because update_level() is used to get a halfstep value in the time evoultion
        if set_flag == True:
            self.set_level(level_new,display_warning=False)
        elif set_flag == False:
            return level_new

    def update_pressure(self,set_flag=False):
        # update pressure based on level and flux velocities 
        # see confluence doc for explaination
        influx_vel  = abs(self.influx/self.area)
        outflux_vel = abs(self.outflux/self.area_out)
        p_new = self.density*self.g*self.level+self.density*outflux_vel*(influx_vel-outflux_vel)
        # set flag for consistency with update_level()
        if set_flag ==True:
            self.set_pressure(p_new,display_warning=False)
        elif set_flag == False:
            return p_new

    def update_volume(self,set_flag=False):
        volume_new = self.level*self.area
        # set flag for consistency with update_level()
        if set_flag == True:
            self.set_volume(volume_new,display_warning=False)
        elif set_flag == False:
            return volume_new

#methods 
    def timestep_reservoir_evolution(self):
        # update outflux, level, pressure and volume based on current pipeline pressure and waterlevel in reservoir
        
        # get some variables
        dt      = self.timestep
        rho     = self.density
        g       = self.g
        A       = self.area
        A_a     = self.area_out
        yn      = self.outflux/A_a # outflux velocity
        h       = self.level
        h_hs    = self.update_level(dt/2)
        p       = self.pressure
        p_hs    = self.pressure + rho*g*(h_hs-h)
        
        # perform explicit 4 step Runge Kutta
        Y1      = yn
        Y2      = yn + dt/2*FODE_function(Y1,h,A,A_a,p,rho,g)
        Y3      = yn + dt/2*FODE_function(Y2,h_hs,A,A_a,p_hs,rho,g)
        Y4      = yn + dt*FODE_function(Y3,h_hs,A,A_a,p_hs,rho,g)
        ynp1    = yn + dt/6*(FODE_function(Y1,h,A,A_a,p,rho,g)+2*FODE_function(Y2,h_hs,A,A_a,p_hs,rho,g)+ \
            2*FODE_function(Y3,h_hs,A,A_a,p_hs,rho,g)+ FODE_function(Y4,h,A,A_a,p,rho,g))

        # set/update the attributes to their new values 
        self.set_outflux(ynp1*A_a,display_warning=False)
        self.update_level(dt,set_flag=True)
        self.update_volume(set_flag=True)
        self.update_pressure(set_flag=True)
