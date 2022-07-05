import numpy as np 

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
    ynp1 = yn + dt/6*(FODE_function(Y1, h, alpha, p0)+2*FODE_function(Y2, h_hs, alpha, p_hs)+ \
        2*FODE_function(Y3, h_hs, alpha, p_hs)+ FODE_function(Y4, h, alpha, p0))

