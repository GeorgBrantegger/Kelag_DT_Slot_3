import numpy as np

#importing pressure conversion function
import sys
import os
current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)
from functions.pressure_conversion import pressure_conversion


class Druckrohrleitung_class:
# units
    acceleration_unit   = r'$\mathrm{m}/\mathrm{s}^2$'
    angle_unit          = '°'
    area_unit           = r'$\mathrm{m}^2$'
    density_unit        = r'$\mathrm{kg}/\mathrm{m}^3$'
    flux_unit           = r'$\mathrm{m}^3/\mathrm{s}$'
    length_unit         = 'm'
    pressure_unit       = 'Pa'
    time_unit           = 's'
    velocity_unit       = r'$\mathrm{m}/\mathrm{s}$' # for flux and pressure propagation
    volume_unit         = r'$\mathrm{m}^3$'

    acceleration_unit_print   = 'm/s²'
    angle_unit_print          = '°'
    area_unit_print           = 'm²'
    density_unit_print        = 'kg/m³'
    flux_unit_print           = 'm³/s'
    length_unit_print         = 'm'
    pressure_unit_print       = 'Pa'
    time_unit_print           = 's'
    velocity_unit_print       = 'm/s' # for flux and pressure propagation
    volume_unit_print         = 'm³'
    

# init
    def __init__(self,total_length,diameter,number_segments,pipeline_angle,Darcy_friction_factor,rho=1000,g=9.81):
        self.length     = total_length
        self.dia        = diameter
        self.n_seg      = number_segments
        self.angle      = pipeline_angle
        self.f_D        = Darcy_friction_factor     # = Rohrreibungszahl oder flow coefficient        
        self.density    = 1000
        self.g          = g
        
        self.dx      = total_length/number_segments
        self.l_vec      = np.arange(0,(number_segments+1)*self.dx,self.dx)

        # initialize for get_info method
        self.c          = '--'
        self.dt         = '--'

# setter
    def set_pressure_propagation_velocity(self,c):
        self.c      = c
        self.dt     = self.dx/c

    def set_number_of_timesteps(self,number_timesteps):
        self.nt     = number_timesteps
        if self.c == '--':
            raise Exception('Please set the pressure propagation velocity before setting the number of timesteps.')
        else:
            self.t_vec = np.arange(0,self.nt*self.dt,self.dt)

    def set_initial_pressure(self,pressure,input_unit = 'Pa'):
        p,_ = pressure_conversion(pressure,input_unit,target_unit=self.pressure_unit)
        if np.size(p) == 1:
            self.p0 = np.full_like(self.l_vec,p)
        elif np.size(p) == np.size(self.l_vec):
            self.p0 = p
        else:
             raise Exception('Unable to assign initial pressure. Input has to be of size 1 or' + np.size(self.l_vec))
        
        #initialize the vectors in which the old and new pressures are stored for the method of characteristics
        self.p_old = self.p0.copy()
        self.p = np.empty_like(self.p_old)
    
    def set_initial_flow_velocity(self,velocity):
        if np.size(velocity) == 1:
            self.v0 = np.full_like(self.l_vec,velocity)
        elif np.size(velocity) == np.size(self.l_vec):
            self.v0 = velocity
        else:
             raise Exception('Unable to assign initial velocity. Input has to be of size 1 or' + np.size(self.l_vec))    

        #initialize the vectors in which the old and new velocities are stored for the method of characteristics
        self.v_old = self.v0.copy()
        self.v = np.empty_like(self.v_old)
    
    def set_boundary_conditions_next_timestep(self,v_reservoir,p_reservoir,v_turbine,input_unit_pressure = 'Pa'):
        rho                 = self.density
        c                   = self.c
        f_D                 = self.f_D
        dt                  = self.dt
        D                   = self.dia
        p_old               = self.p_old[-2]    # @ second to last node (the one before the turbine)
        v_old               = self.v_old[-2]    # @ second to last node (the one before the turbine)
        self.v_boundary_res = v_reservoir
        self.v_boundary_tur = v_turbine
        self.p_boundary_res,_ = pressure_conversion(p_reservoir,input_unit_pressure,target_unit=self.pressure_unit)
        self.p_boundary_tur = p_old+rho*c*v_old-rho*c*f_D*dt/(2*D)*abs(v_old)*v_old
        self.v[0]       = self.v_boundary_res.copy()
        self.v[-1]      = self.v_boundary_tur.copy()
        self.p[0]       = self.p_boundary_res.copy()
        self.p[-1]      = self.p_boundary_tur.copy()

# getter
    def get_info(self):
        new_line = '\n'

        # :<10 pads the self.value to be 10 characters wide
        print_str = (f"The pipeline has the following attributes: {new_line}" 
            f"----------------------------- {new_line}"
            f"Length                =       {self.length:<10} {self.length_unit_print} {new_line}"
            f"Diameter              =       {self.dia:<10} {self.length_unit_print} {new_line}"
            f"Number of segemnts    =       {self.n_seg:<10} {new_line}"
            f"Number of nodes       =       {self.n_seg+1:<10} {new_line}"
            f"Length per segment    =       {self.dx:<10} {self.length_unit_print} {new_line}"
            f"Pipeline angle        =       {self.angle:<10} {self.angle_unit_print} {new_line}"
            f"Darcy friction factor =       {self.f_D:<10} {new_line}"
            f"Density of liquid     =       {self.density:<10} {self.density_unit_print} {new_line}"
            f"Pressure wave vel.    =       {self.c:<10} {self.velocity_unit_print} {new_line}"
            f"Simulation timesteps  =       {self.dt:<10} {self.time_unit_print } {new_line}"
            f"Number of timesteps   =       {self.nt:<10} {new_line}"
            f"----------------------------- {new_line}"
            f"Velocity and pressure distribution are vectors and are accessible by the .v and .p attribute of the pipeline object")

        print(print_str)    
        

    def get_boundary_conditions_next_timestep(self,target_unit_pressure ='bar'):
        print('The pressure at the reservoir for the next timestep is', '\n', \
                pressure_conversion(self.p_boundary_res,self.pressure_unit_print,target_unit_pressure), '\n', \
            'The velocity at the reservoir for the next timestep is', '\n', \
                self.v_boundary_res, self.velocity_unit, '\n', \
            'The pressure at the turbine for the next timestep is', '\n', \
                pressure_conversion(self.p_boundary_tur,self.pressure_unit_print,target_unit_pressure), '\n', \
            'The velocity at the turbine for the next timestep is', '\n', \
                self.v_boundary_tur, self.velocity_unit)         


    def timestep_characteristic_method(self):
        #number of nodes
        nn  = self.n_seg+1
        rho = self.density
        c   = self.c
        f_D = self.f_D
        dt  = self.dt
        D   = self.dia

        for i in range(1,nn-1):
            self.v[i] = 0.5*(self.v_old[i-1]+self.v_old[i+1])+0.5/(rho*c)*(self.p_old[i-1]-self.p_old[i+1]) \
                -f_D*dt/(4*D)*(abs(self.v_old[i-1])*self.v_old[i-1]+abs(self.v_old[i+1])*self.v_old[i+1])

            self.p[i] = 0.5*rho*c*(self.v_old[i-1]-self.v_old[i+1])+0.5*(self.p_old[i-1]+self.p_old[i+1]) \
                -rho*c*f_D*dt/(4*D)*(abs(self.v_old[i-1])*self.v_old[i-1]-abs(self.v_old[i+1])*self.v_old[i+1])

        self.p_old = self.p.copy()
        self.v_old = self.v.copy()        











