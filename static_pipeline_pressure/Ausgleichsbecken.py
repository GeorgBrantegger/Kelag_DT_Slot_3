import numpy as np 

def Volume_trend(influx, outflux, timestep=1, V_0=0):
    '''
    Returns the trend and the volume and the final volume,  defined
    by influx and outflux patterns. The optional parameter timestep
    defines the time increment over which the fluxes are changing.
    '''
    net_flux = influx-outflux
    delta_V = net_flux*timestep
    V_trend = V_0+np.cumsum(delta_V)
    V_end = V_trend[-1]
    return V_end,  V_trend

def Height_trend(V_trend, area=1, h_crit_low=-np.inf, h_crit_high=np.inf):
    '''
    Returns the trend and the height and the final height,  defined
    by influx and outflux patterns as well as the crosssection area.
    The optional parameters h_crit_low/high indicate limits that the height
    should never exceed. If this occures,  TRUE is returned in the corresponding
    h_crit_flag.
    '''
    h_trend = V_trend/area
    h_crit_flag_low = np.any(h_trend <= h_crit_low)
    h_crit_flag_high = np.any(h_trend >= h_crit_high)
    h_end = h_trend[-1]
    return h_trend, h_end, h_crit_flag_low, h_crit_flag_high

def get_h_halfstep(initial_height, influx, outflux, timestep, area):
    h0      = initial_height
    Q_in    = influx
    Q_out   = outflux
    dt      = timestep
    A       = area

    h_halfstep = h0+1/A*(Q_in-Q_out)*dt/2

def get_p_halfstep(p0, p1):
    p_halfstep = (p0+p1)/2

def FODE_function(x, h, alpha, p, rho=1000., g=9.81):
    f = x*abs(x)/h*alpha+g-p/(rho*h)
    return f


def e_RK_4(yn, h, dt, Q0, Q1, A0, A1, p0, p1):
    alpha = (A1/A0-1)

    h_hs = get_h_halfstep(h, Q0, Q1, dt, A0)
    p_hs = get_p_halfstep(p0, p1)
    Y1 = yn
    Y2 = yn + dt/2*FODE_function(Y1, h, alpha, p0)
    Y3 = yn + dt/2*FODE_function(Y2, h_hs, alpha, p_hs)
    Y4 = yn + dt*FODE_function(Y3, h_hs, alpha, p_hs)
    ynp1 = yn + dt/6*(FODE_function(Y1, h, alpha, p)+2*FODE_function(Y2, h_hs, alpha, p_hs)+ \
        2*FODE_function(Y3, h_hs, alpha, p_hs)+ FODE_function(Y4, h, alpha, p))




## testing
# if __name__ == "__main__":
#     influx = np.full([1, 100],  6)
#     outflux = np.full_like(influx,  4)
#     V_end,  V_trend = Volume_trend(influx,  outflux, timestep=0.5, V_0 = 100)
#     print(V_end)
#     print(V_trend)
