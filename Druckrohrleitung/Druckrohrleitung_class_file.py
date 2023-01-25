import os
import sys

import numpy as np

#importing pressure conversion function
current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)
from functions.pressure_conversion import pressure_conversion


class Druckrohrleitung_class:
# units 
    # make sure that units and display units are the same
    # units are used to label graphs and disp units are used to have a bearable format when using pythons print()
    acceleration_unit           = r'$\mathrm{m}/\mathrm{s}^2$'
    angle_unit                  = 'rad'
    area_unit                   = r'$\mathrm{m}^2$'
    density_unit                = r'$\mathrm{kg}/\mathrm{m}^3$'
    flux_unit                   = r'$\mathrm{m}^3/\mathrm{s}$'
    length_unit                 = 'm'
    pressure_unit               = 'Pa'  # DONT CHANGE needed for pressure conversion
    time_unit                   = 's'
    velocity_unit               = r'$\mathrm{m}/\mathrm{s}$' # for flux and pressure propagation
    volume_unit                 = r'$\mathrm{m}^3$'

    acceleration_unit_disp      = 'm/s²'
    angle_unit_disp             = 'rad'
    area_unit_disp              = 'm²'
    density_unit_disp           = 'kg/m³'
    flux_unit_disp              = 'm³/s'
    length_unit_disp            = 'm'
    # pressure_unit_disp will be set within the __init__() method
    time_unit_disp              = 's'
    velocity_unit_disp          = 'm/s' # for flux and pressure propagation
    volume_unit_disp            = 'm³'

    g = 9.81    # m/s²  gravitational acceleration

# init
    def __init__(self,total_length,diameter,pipeline_head,number_segments,Darcy_friction_factor,pw_vel,timestep,pressure_unit_disp,rho=1000):
        """
        Creates a reservoir with given attributes in this order: \n
            Pipeline length                 [m]         \n
            Pipeline diameter               [m]         \n
            Pipeline head                   [m]         \n
            Number of pipeline segments     [1]         \n
            Darcy friction factor           [1]         \n
            Pressure wave velocity          [m/s]       \n
            Simulation timestep             [s]         \n
            Pressure unit for displaying    [string]    \n            
            Density of the liquid           [kg/m³]     \n
        """
        self.length     = total_length  	                                        # total length of the pipeline
        self.dia        = diameter                                                  # diameter of the pipeline
        self.head       = pipeline_head                                             # hydraulic head of the pipeline
        self.n_seg      = number_segments                                           # number of segments for the method of characteristics
        self.f_D        = Darcy_friction_factor                                     # = Rohrreibungszahl oder flow coefficient  
        self.c          = pw_vel                                                    # propagation velocity of pressure wave
        self.dt         = timestep
        self.density    = rho                                                       # density of the liquid in the pipeline
        
        # derivatives of input attributes
        self.angle      = np.arcsin(self.head/self.length)                          # angle of the pipeline
        self.A          = (diameter/2)**2*np.pi                                     # crossectional area of the pipeline
        self.dx         = total_length/number_segments                              # length of each segment
        self.x_vec      = np.arange(0,(number_segments+1),1)*self.dx                # vector giving the distance from each node to the start of the pipeline
        self.h_vec      = np.arange(0,(number_segments+1),1)*self.head/self.n_seg   # vector giving the height difference from each node to the start of the pipeline
        
        self.pressure_unit_disp = pressure_unit_disp                                # pressure unit for displaying

# setter - set attributes
    def set_initial_pressure(self,pressure,display_warning=True):
        # initialize the pressure distribution in the pipeline
        if display_warning == True:
            print('You are setting the pressure distribution in the pipeline manually. \n \
                This is not an intended use of this method. \n \
                Refer to the set_steady_state() method instead.')

        # make sure the vector has the right size
        if np.size(pressure)   == 1:
            p0 = np.full_like(self.x_vec,pressure)
        elif np.size(pressure) == np.size(self.x_vec):
            p0 = pressure
        else:
             raise Exception('Unable to assign initial pressure. Input has to be of size 1 or' + np.size(self.x_vec))
        
        #initialize the vectors in which the old and new pressures are stored for the method of characteristics
        self.p_old  = p0.copy()
        self.p      = p0.copy()
        self.p0     = p0.copy()
        # initialize the vectors in which the minimal and maximal value of the pressure at each node are stores
        self.p_min  = p0.copy()
        self.p_max  = p0.copy()
    
    def set_initial_flow_velocity(self,velocity, display_warning=True):
        # initialize the velocity distribution in the pipeline
        if display_warning == True:
            print('You are setting the velocity distribution in the pipeline manually. \n \
                This is not an intended use of this method. \n \
                Refer to the set_steady_state() method instead.')   

        # make sure the vector has the right size  
        if np.size(velocity)   == 1:
            v0 = np.full_like(self.x_vec,velocity)
        elif np.size(velocity) == np.size(self.x_vec):
            v0 = velocity
        else:
             raise Exception('Unable to assign initial velocity. Input has to be of size 1 or' + np.size(self.x_vec))    

        #initialize the vectors in which the old and new velocities are stored for the method of characteristics
        self.v_old  = v0.copy()
        self.v      = v0.copy()
        # initialize the vectors in which the minimal and maximal value of the velocity at each node are stores
        self.v_min  = v0.copy()
        self.v_max  = v0.copy()
    
    def set_boundary_conditions_next_timestep(self,p_reservoir,v_turbine):
        # derived from the method of characteristics, one can set the boundary conditions for the pressures and flow velocities at the reservoir and the turbine
        # the boundary velocity at the turbine is specified by the flux through the turbine or an external boundary condition
            # the pressure at the turbine will be calculated using the forward characteristic
        # the boundary pressure at the reservoir is specified by the level in the reservoir of an external boundary condition
            # the velocity at the reservoir will be calculated using the backward characteristic

        # constants for a cleaner formula
        rho                 = self.density
        c                   = self.c
        f_D                 = self.f_D
        dt                  = self.dt
        D                   = self.dia
        g                   = self.g
        alpha               = self.angle
        p_old_tur           = self.p_old[-2]    # @ second to last node (the one before the turbine)
        v_old_tur           = self.v_old[-2]    # @ second to last node (the one before the turbine)
        p_old_res           = self.p_old[1]     # @ second  node (the one after the reservoir)
        v_old_res           = self.v_old[1]     # @ second  node (the one after the reservoir)
        # set the boundary conditions derived from reservoir and turbine
        v_boundary_tur      = v_turbine         # at new timestep
        p_boundary_res      = p_reservoir       # at new timestep
        # calculate the missing boundary conditions
        v_boundary_res      = v_old_res+1/(rho*c)*(p_boundary_res-p_old_res)+dt*g*np.sin(alpha)-f_D*dt/(2*D)*abs(v_old_res)*v_old_res
        p_boundary_tur      = p_old_tur-rho*c*(v_boundary_tur-v_old_tur)+rho*c*dt*g*np.sin(alpha)-f_D*rho*c*dt/(2*D)*abs(v_old_tur)*v_old_tur

        # write boundary conditions to the velocity/pressure vectors of the next timestep
        self.v[0]           = v_boundary_res
        self.v[-1]          = v_boundary_tur
        self.p[0]           = p_boundary_res
        self.p[-1]          = p_boundary_tur

    def set_steady_state(self,ss_flux,ss_pressure_res):
        # set the pressure and velocity distributions, that allow a constant flow of water from the (steady-state) reservoir to the (steady-state) turbine
            # the flow velocity is given by the constant flow through the pipe
        ss_v0 = np.full_like(self.x_vec,ss_flux/self.A)

        # the static pressure is given by static state pressure of the reservoir, corrected for the hydraulic head of the pipe and friction losses
        ss_pressure     = ss_pressure_res+(self.density*self.g*self.h_vec)-(self.f_D*self.x_vec/self.dia*self.density/2*ss_v0**2)

        # set the initial conditions
        self.set_initial_flow_velocity(ss_v0,display_warning=False)
        self.set_initial_pressure(ss_pressure,display_warning=False)

# getter - return attributes
    def get_info(self):
        new_line    = '\n'
        angle_deg   = round(self.angle/np.pi*180,3)


        # :<10 pads the self.value to be 10 characters wide
        print_str = (f"The pipeline has the following attributes: {new_line}" 
            f"----------------------------- {new_line}"
            f"Length                =       {self.length:<10} {self.length_unit_disp} {new_line}"
            f"Diameter              =       {self.dia:<10} {self.length_unit_disp} {new_line}"
            f"Hydraulic head        =       {self.head:<10} {self.length_unit_disp} {new_line}"
            f"Number of segments    =       {self.n_seg:<10} {new_line}"
            f"Number of nodes       =       {self.n_seg+1:<10} {new_line}"
            f"Length per segments   =       {self.dx:<10} {self.length_unit_disp} {new_line}"
            f"Pipeline angle        =       {round(self.angle,3):<10} {self.angle_unit_disp} {new_line}"
            f"Pipeline angle        =       {angle_deg}° {new_line}"
            f"Darcy friction factor =       {self.f_D:<10} {new_line}"
            f"Density of liquid     =       {self.density:<10} {self.density_unit_disp} {new_line}"
            f"Pressure wave vel.    =       {self.c:<10} {self.velocity_unit_disp} {new_line}"
            f"Simulation timestep   =       {self.dt:<10} {self.time_unit_disp} {new_line}"
            f"----------------------------- {new_line}"
            f"Velocity and pressure distribution are vectors and are accessible by the .v and .p attribute of the pipeline object")

        print(print_str)    

    def get_current_pressure_distribution(self,disp_flag=False):
        # disp_flag if one wants to directly plot the return of this method
        if disp_flag == True:       # convert to pressure unit disp
            return pressure_conversion(self.p,self.pressure_unit,self.pressure_unit_disp)
        elif disp_flag == False:    # stay in Pa
            return self.p

    def get_current_velocity_distribution(self):
        return self.v

    def get_current_flux_distribution(self):
        return self.v*self.A

    def get_lowest_pressure_per_node(self,disp_flag=False):
        if disp_flag == True:       # convert to pressure unit disp
            return pressure_conversion(self.p_min,self.pressure_unit,self.pressure_unit_disp)
        elif disp_flag == False:    # stay in Pa
            return self.p_min

    def get_highest_pressure_per_node(self,disp_flag=False):
        if disp_flag == True:       # convert to pressure unit disp
            return pressure_conversion(self.p_max,self.pressure_unit,self.pressure_unit_disp)
        elif disp_flag == False:    # stay in Pa
            return self.p_max

    def get_lowest_velocity_per_node(self):
        return self.v_min

    def get_highest_velocity_per_node(self):
        return self.v_max

    def get_lowest_flux_per_node(self):
        return self.v_min*self.A

    def get_highest_flux_per_node(self):
        return self.v_max*self.A

    def get_initial_pressure_distribution(self,disp_flag=False):
        # disp_flag if one wants to directly plot the return of this method
        if disp_flag == True:       # convert to pressure unit disp
            return pressure_conversion(self.p0,self.pressure_unit,self.pressure_unit_disp)
        elif disp_flag == False:    # stay in Pa
            return self.p0    

    def timestep_characteristic_method(self):
    # use the method of characteristics to calculate the pressure and velocities at all nodes except the boundary ones
        # they are set with the .set_boundary_conditions_next_timestep() method beforehand
        
        # constants for cleaner formula
        nn      = self.n_seg+1      # number of nodes
        rho     = self.density      # density of liquid
        c       = self.c            # pressure propagation velocity
        f_D     = self.f_D          # Darcy friction coefficient
        dt      = self.dt           # timestep
        D       = self.dia          # pipeline diameter
        g       = self.g            # graviational acceleration
        alpha   = self.angle        # pipeline angle

        # Vectorize this loop?
        for i in range(1,nn-1):
            self.v[i] = 0.5*(self.v_old[i+1]+self.v_old[i-1])-0.5/(rho*c)*(self.p_old[i+1]-self.p_old[i-1]) \
                +dt*g*np.sin(alpha)-f_D*dt/(4*D)*(abs(self.v_old[i+1])*self.v_old[i+1]+abs(self.v_old[i-1])*self.v_old[i-1])

            self.p[i] = 0.5*(self.p_old[i+1]+self.p_old[i-1])-0.5*rho*c*(self.v_old[i+1]-self.v_old[i-1]) \
                +f_D*rho*c*dt/(4*D)*(abs(self.v_old[i+1])*self.v_old[i+1]-abs(self.v_old[i-1])*self.v_old[i-1])

        # update overall min and max values for pressure and velocity per node
        self.p_min = np.minimum(self.p_min,self.p)
        self.p_max = np.maximum(self.p_max,self.p)
        self.v_min = np.minimum(self.v_min,self.v)
        self.v_max = np.maximum(self.v_max,self.v)

        # prepare for next call
            # use .copy() to write data to another memory location and avoid the usual python reference pointer
            # else one can overwrite data by accidient and change two variables at once without noticing
        self.p_old = self.p.copy()
        self.v_old = self.v.copy()        

    def timestep_characteristic_method_vectorized(self):
    # use the method of characteristics to calculate the pressure and velocities at all nodes except the boundary ones
        # they are set with the .set_boundary_conditions_next_timestep() method beforehand
        
        # constants for cleaner formula
        rho     = self.density      # density of liquid
        c       = self.c            # pressure propagation velocity
        f_D     = self.f_D          # Darcy friction coefficient
        dt      = self.dt           # timestep
        D       = self.dia          # pipeline diameter
        g       = self.g            # graviational acceleration
        alpha   = self.angle        # pipeline angle

        # Vectorized loop
        self.v[1:-1] = 0.5*(self.v_old[2:]+self.v_old[:-2])-0.5/(rho*c)*(self.p_old[2:]-self.p_old[:-2]) \
            +dt*g*np.sin(alpha)-f_D*dt/(4*D)*(np.abs(self.v_old[2:])*self.v_old[2:]+np.abs(self.v_old[:-2])*self.v_old[:-2])

        self.p[1:-1] = 0.5*(self.p_old[2:]+self.p_old[:-2])-0.5*rho*c*(self.v_old[2:]-self.v_old[:-2]) \
            +f_D*rho*c*dt/(4*D)*(np.abs(self.v_old[2:])*self.v_old[2:]-np.abs(self.v_old[:-2])*self.v_old[:-2])

        # update overall min and max values for pressure and velocity per node
        self.p_min = np.minimum(self.p_min,self.p)
        self.p_max = np.maximum(self.p_max,self.p)
        self.v_min = np.minimum(self.v_min,self.v)
        self.v_max = np.maximum(self.v_max,self.v)

        # prepare for next call
            # use .copy() to write data to another memory location and avoid the usual python reference pointer
            # else one can overwrite data by accidient and change two variables at once without noticing
        self.p_old = self.p.copy()
        self.v_old = self.v.copy()        
