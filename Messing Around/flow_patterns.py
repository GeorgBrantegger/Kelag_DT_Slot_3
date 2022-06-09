import numpy as np
import pandas as pd
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go


def return_flux_profiles(number_of_steps = 1,influx_identifier = 0, outflux_identifier = 0,influx_offset=0,outflux_offset=0, outflux_delay = 0):
    ''' Identifier patterns:
    0                   ... constant
    'lin_SSSS'          ... linear increase with slope int(SSSS)
    'st_SSSS_PPPP'      ... sawtooth pattern with slope int(SSSS) and period int(PPPP) steps
    '''
    
    # case identifiers for if statment
    i = influx_identifier
    o = outflux_identifier
    
    n = number_of_steps
    #starting value for the influx and outflux
    i_o = influx_offset
    o_o = outflux_offset
    # number of steps, the outflux is held at 0 at the beginning
    o_d = outflux_delay


# get base profile for the influx (offset will get applied later)
    if i == 0:
        influx_profile = np.zeros(n)
    elif 'lin' in influx_identifier:
        k = int(influx_identifier[-4:])
        influx_profile = np.linspace(0,k*(n-1),n)
    elif 'st' in influx_identifier:
        k = int(influx_identifier[3:7])
        p = int(influx_identifier[-4:])
        influx_profile = np.tile(np.linspace(0,k*(p-1),p),int(np.ceil(n/p)))


# apply influx offset
    influx_profile = influx_offset + influx_profile

    if o == 0:
        outflux_profile = np.zeros(n)
    elif 'lin' in outflux_identifier:
        k = int(outflux_identifier[-4:])
        outflux_profile = np.linspace(0,k*(n-1),n)
    elif 'st' in outflux_identifier:
        k = int(outflux_identifier[3:7])
        p = int(outflux_identifier[-4:])
        outflux_profile = np.tile(np.linspace(0,k*(p-1),p),int(np.ceil(n/p)))

#apply outflux offset and delay (delay means, that the first o_d steps, the outflux will be 0)
    outflux_profile = np.concatenate((np.zeros(o_d),outflux_profile[:-o_d]+o_o))
    
    return influx_profile,outflux_profile

def make_flux_df(influx_profile,outflux_profile, time = 0):
    if time == 0:
        time = np.arange(0,len(influx_profile))
    flux_df = pd.DataFrame(np.transpose([time, influx_profile, outflux_profile]), \
        columns=['time', 'influx', 'outflux'])
    return flux_df


if __name__ == "__main__":
    influx_profile,outflux_profile = return_flux_profiles(100,influx_identifier='st_0010_0010',influx_offset=10)
    print(influx_profile)