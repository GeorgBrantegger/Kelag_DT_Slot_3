from pressure_conversion import pressure_conversion
import numpy as np

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

        # workaround for try-except construct in set_number_of_timesteps
        self.c          = 0

# setter
    def set_pressure_propagation_velocity(self,c):
        self.c      = c
        self.dt     = self.dx/c

    def set_number_of_timesteps(self,number_timesteps):
        self.nt     = number_timesteps
        if self.c == 0:
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
        self.p_new = np.empty_like(self.p_old)
    
    def set_initial_flow_velocity(self,velocity):
        if np.size(velocity) == 1:
            self.v0 = np.full_like(self.l_vec,velocity)
        elif np.size(velocity) == np.size(self.l_vec):
            self.v0 = velocity
        else:
             raise Exception('Unable to assign initial velocity. Input has to be of size 1 or' + np.size(self.l_vec))    

        #initialize the vectors in which the old and new velocities are stored for the method of characteristics
        self.v_old = self.v0.copy()
        self.v_new = np.empty_like(self.v_old)
    
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
        self.v_new[0]       = self.v_boundary_res.copy()
        self.v_new[-1]      = self.v_boundary_tur.copy()
        self.p_new[0]       = self.p_boundary_res.copy()
        self.p_new[-1]      = self.p_boundary_tur.copy()

# getter
    def get_pipeline_geometry(self):
        print('The total length of the pipeline is', '\n', \
             self.length, self.length_unit, '\n', \
            'The diameter of the pipeline is', '\n', \
            self.dia, self.length_unit, '\n', \
            'The pipeline is divided into', self.n_seg , 'segments of length', '\n', \
            round(self.dx,1), self.length_unit, '\n', \
            'The pipeline has an inclination angle of', '\n', \
            self.angle, self.angle_unit)
    
    def get_other_pipeline_info(self):
        print('The Darcy-friction factor of the pipeline is', '\n', \
            self.f_D, '\n', \
            'The pipeline is filled with a liquid with density', '\n', \
            self.density, self.density_unit, '\n', \
            'The gravitational acceleration is set to', '\n', \
            self.g, self.acceleration_unit)
    
    def get_pressure_propagation_velocity(self):
        print('The pressure propagation velocity in the pipeline is', '\n', \
            self.c, self.velocity_unit)

    def get_number_of_timesteps(self):
        print(self.nt, 'timesteps are performed in the simulation')

    
    def get_initial_pressure(self,target_unit='bar'):
        print('The inital pressure distribution in is', '\n', \
            pressure_conversion(self.p0,self.pressure_unit,target_unit))

    def get_initial_flow_velocity(self):  
        print('The inital velocity distribution is', '\n', \
            self.v0, self.velocity_unit)        

    def get_boundary_conditions_next_timestep(self,target_unit_pressure ='bar'):
        print('The pressure at the reservoir for the next timestep is', '\n', \
                pressure_conversion(self.p_boundary_res,self.pressure_unit,target_unit_pressure), '\n', \
            'The velocity at the reservoir for the next timestep is', '\n', \
                self.v_boundary_res, self.velocity_unit, '\n', \
            'The pressure at the turbine for the next timestep is', '\n', \
                pressure_conversion(self.p_boundary_tur,self.pressure_unit,target_unit_pressure), '\n', \
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
            self.v_new[i] = 0.5*(self.v_old[i-1]+self.v_old[i+1])+0.5/(rho*c)*(self.p_old[i-1]-self.p_old[i+1]) \
                -f_D*dt/(4*D)*(abs(self.v_old[i-1])*self.v_old[i-1]+abs(self.v_old[i+1])*self.v_old[i+1])

            self.p_new[i] = 0.5*rho*c*(self.v_old[i-1]-self.v_old[i+1])+0.5*(self.p_old[i-1]+self.p_old[i+1]) \
                -rho*c*f_D*dt/(4*D)*(abs(self.v_old[i-1])*self.v_old[i-1]-abs(self.v_old[i+1])*self.v_old[i+1])

        self.p_old = self.p_new.copy()
        self.v_old = self.v_new.copy()        











