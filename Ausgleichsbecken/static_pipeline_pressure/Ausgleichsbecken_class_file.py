from Ausgleichsbecken import FODE_function, get_h_halfstep, get_p_halfstep
from pressure_conversion import pressure_conversion
class Ausgleichsbecken_class:
# units
    area_unit           = r'$\mathrm{m}^2$'
    area_outflux_unit   = r'$\mathrm{m}^2$'
    level_unit          = 'm'
    volume_unit         = r'$\mathrm{m}^3$'
    flux_unit           = r'$\mathrm{m}^3/\mathrm{s}$'
    time_unit           = 's'
    pressure_unit       = 'Pa'

# init
    def __init__(self,area,outflux_area,level_min,level_max,timestep = 1):
        self.area           = area              # base area of the rectangular structure
        self.area_outflux   = outflux_area      # area of the outlet towards the pipeline        
        self.level_min      = level_min         # lowest  allowed water level    
        self.level_max      = level_max         # highest allowed water level
        self.timestep       = timestep          # timestep of the simulation    

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
    def get_area(self):
        print('The base area of the cuboid reservoir is', self.area, self.area_unit)

    def get_outflux_area(self):
        print('The outflux area from the cuboid reservoir to the pipeline is', \
            self.area_outflux, self.area_outflux_unit)
    
    def get_level(self):
        print('The current level in the reservoir is', self.level , self.level_unit)

    def get_crit_levels(self):
        print('The critical water levels in the reservoir are:  \n',\
            '   Minimum:', self.level_min , self.level_unit , '\n',\
            '   Maximum:', self.level_max , self.level_unit )

    def get_volume(self):
        print('The current water volume in the reservoir is', self.volume, self.volume_unit)

    def get_timestep(self):
        print('The timestep for the simulation is' , self.timestep, self.time_unit)

    def get_influx(self):
        print('The current influx is', self.influx, self.flux_unit)

    def get_outflux(self):
        print('The current outflux is', self.outflux, self.flux_unit)

# methods
    def update_level(self,timestep):
        net_flux = self.influx-self.outflux
        delta_V = net_flux*timestep
        new_level = (self.volume+delta_V)/self.area
        return new_level


    def e_RK_4(self):
        # Update to deal with non constant pipeline pressure!
        yn = self.outflux/self.area_outflux
        h = self.level
        dt = self.timestep
        p,_ = pressure_conversion(self.initial_pressure,self.pressure_unit,'Pa')
        p_hs,_ = pressure_conversion(self.initial_pressure,self.pressure_unit,'Pa')
        alpha = (self.area_outflux/self.area-1)
        h_hs = self.update_level(dt/2)
        Y1 = yn
        Y2 = yn + dt/2*FODE_function(Y1, h, alpha, self.initial_pressure)
        Y3 = yn + dt/2*FODE_function(Y2, h_hs, alpha, p_hs)
        Y4 = yn + dt*FODE_function(Y3, h_hs, alpha, p_hs)
        ynp1 = yn + dt/6*(FODE_function(Y1, h, alpha, p)+2*FODE_function(Y2, h_hs, alpha, p_hs)+ \
            2*FODE_function(Y3, h_hs, alpha, p_hs)+ FODE_function(Y4, h, alpha, p))

        self.outflux = ynp1*self.area_outflux