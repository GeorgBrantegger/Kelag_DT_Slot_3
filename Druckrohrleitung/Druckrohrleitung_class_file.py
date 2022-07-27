import numpy as np

class Druckrohrleitung_class:
# units
    acceleration_unit           = r'$\mathrm{m}/\mathrm{s}^2$'
    angle_unit                  = 'rad'
    area_unit                   = r'$\mathrm{m}^2$'
    density_unit                = r'$\mathrm{kg}/\mathrm{m}^3$'
    flux_unit                   = r'$\mathrm{m}^3/\mathrm{s}$'
    length_unit                 = 'm'
    pressure_unit               = 'Pa'
    time_unit                   = 's'
    velocity_unit               = r'$\mathrm{m}/\mathrm{s}$' # for flux and pressure propagation
    volume_unit                 = r'$\mathrm{m}^3$'

    acceleration_unit_print     = 'm/s²'
    angle_unit_print            = 'rad'
    area_unit_print             = 'm²'
    density_unit_print          = 'kg/m³'
    flux_unit_print             = 'm³/s'
    length_unit_print           = 'm'
    time_unit_print             = 's'
    velocity_unit_print         = 'm/s' # for flux and pressure propagation
    volume_unit_print           = 'm³'

# init
    def __init__(self,total_length,diameter,number_segments,pipeline_angle,Darcy_friction_factor,rho=1000,g=9.81):
        self.length     = total_length  	                            # total length of the pipeline
        self.dia        = diameter                                      # diameter of the pipeline
        self.n_seg      = number_segments                               # number of segments for the method of characteristics
        self.angle      = pipeline_angle                                # angle of the pipeline
        self.f_D        = Darcy_friction_factor                         # = Rohrreibungszahl oder flow coefficient        
        self.density    = rho                                           # density of the liquid in the pipeline
        self.g          = g                                             # gravitational acceleration
        
        self.A          = (diameter/2)**2*np.pi

        self.dx         = total_length/number_segments                  # length of each segment
        self.l_vec      = np.arange(0,(number_segments+1),1)*self.dx    # vector giving the distance from each node to the start of the pipeline

        # initialize for get_info method
        self.c          = '--'
        self.dt         = '--'

# setter
    def set_pressure_propagation_velocity(self,c):
        self.c      = c         # propagation velocity of the pressure wave
        self.dt     = self.dx/c # timestep derived from c, demanded by the method of characteristics

    def set_number_of_timesteps(self,number_timesteps):
        self.nt     = number_timesteps  # number of timesteps 
        if self.c == '--':
            raise Exception('Please set the pressure propagation velocity before setting the number of timesteps.')
        else:
            self.t_vec = np.arange(0,self.nt*self.dt,self.dt)

    def set_initial_pressure(self,pressure):
        # initialize the pressure distribution in the pipeline
        if np.size(pressure)   == 1:
            self.p0 = np.full_like(self.l_vec,pressure)
        elif np.size(pressure) == np.size(self.l_vec):
            self.p0 = pressure
        else:
             raise Exception('Unable to assign initial pressure. Input has to be of size 1 or' + np.size(self.l_vec))
        
        #initialize the vectors in which the old and new pressures are stored for the method of characteristics
        self.p_old  = self.p0.copy()
        self.p      = self.p0.copy()
    
    def set_initial_flow_velocity(self,velocity):
        # initialize the velocity distribution in the pipeline
        if np.size(velocity)   == 1:
            self.v0 = np.full_like(self.l_vec,velocity)
        elif np.size(velocity) == np.size(self.l_vec):
            self.v0 = velocity
        else:
             raise Exception('Unable to assign initial velocity. Input has to be of size 1 or' + np.size(self.l_vec))    

        #initialize the vectors in which the old and new velocities are stored for the method of characteristics
        self.v_old  = self.v0.copy()
        self.v      = self.v0.copy()
    
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

    def set_steady_state(self,ss_flux,ss_level_reservoir,pl_vec,h_vec):
        # set the pressure and velocity distributions, that allow a constant flow of water from the (steady-state) reservoir to the (steady-state) turbine
            # the flow velocity is given by the constant flow through the pipe
        ss_v0 = np.full(self.n_seg+1,ss_flux/self.A)
            # the static pressure is given by the hydrostatic pressure, corrected for friction losses and dynamic pressure
        ss_pressure = self.density*self.g*(ss_level_reservoir+h_vec)-ss_v0**2*self.density/2-(self.f_D*pl_vec/self.dia*self.density/2*ss_v0**2)

        self.set_initial_flow_velocity(ss_v0)
        self.set_initial_pressure(ss_pressure)

# getter
    def get_info(self):
        new_line    = '\n'
        angle_deg   = round(self.angle/np.pi*180,3)


        # :<10 pads the self.value to be 10 characters wide
        print_str = (f"The pipeline has the following attributes: {new_line}" 
            f"----------------------------- {new_line}"
            f"Length                =       {self.length:<10} {self.length_unit_print} {new_line}"
            f"Diameter              =       {self.dia:<10} {self.length_unit_print} {new_line}"
            f"Number of segments    =       {self.n_seg:<10} {new_line}"
            f"Number of nodes       =       {self.n_seg+1:<10} {new_line}"
            f"Length per segments   =       {self.dx:<10} {self.length_unit_print} {new_line}"
            f"Pipeline angle        =       {round(self.angle,3):<10} {self.angle_unit_print} {new_line}"
            f"Pipeline angle        =       {angle_deg}° {new_line}"
            f"Darcy friction factor =       {self.f_D:<10} {new_line}"
            f"Density of liquid     =       {self.density:<10} {self.density_unit_print} {new_line}"
            f"Pressure wave vel.    =       {self.c:<10} {self.velocity_unit_print} {new_line}"
            f"Simulation timestep   =       {self.dt:<10} {self.time_unit_print} {new_line}"
            f"Number of timesteps   =       {self.nt:<10} {new_line}"
            f"Total simulation time =       {self.nt*self.dt:<10} {self.time_unit_print} {new_line}"
            f"----------------------------- {new_line}"
            f"Velocity and pressure distribution are vectors and are accessible by the .v and .p attribute of the pipeline object")

        print(print_str)    

    def get_current_pressure_distribution(self):
        return self.p

    def get_current_velocity_distribution(self):
        return self.v

    def timestep_characteristic_method(self):
    # use the method of characteristics to calculate the pressure and velocities at all nodes except the boundary ones
        # they are set with the  .set_boundary_conditions_next_timestep() method beforehand
        
        nn      = self.n_seg+1      # number of nodes
        rho     = self.density      # density of liquid
        c       = self.c            # pressure propagation velocity
        f_D     = self.f_D          # Darcy friction coefficient
        dt      = self.dt           # timestep
        D       = self.dia          # pipeline diametet
        g       = self.g            # graviational acceleration
        alpha   = self.angle        # pipeline angle

        # Vectorize this loop?
        for i in range(1,nn-1):
            self.v[i] = 0.5*(self.v_old[i+1]+self.v_old[i-1])-0.5/(rho*c)*(self.p_old[i+1]-self.p_old[i-1]) \
                +dt*g*np.sin(alpha)-f_D*dt/(4*D)*(abs(self.v_old[i+1])*self.v_old[i+1]+abs(self.v_old[i-1])*self.v_old[i-1])

            self.p[i] = 0.5*(self.p_old[i+1]+self.p_old[i-1]) - 0.5*rho*c*(self.v_old[i+1]-self.v_old[i-1]) \
                +f_D*rho*c*dt/(4*D)*(abs(self.v_old[i+1])*self.v_old[i+1]-abs(self.v_old[i-1])*self.v_old[i-1])

        # prepare for next call
            # use .copy() to write data to another memory location and avoid the usual python reference pointer
            # else one can overwrite data by accidient and change two variables at once without noticing
        self.p_old = self.p.copy()
        self.v_old = self.v.copy()        
