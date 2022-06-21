import numpy as np 

def simple_time_delay(delta_p_profile,delay,timestep):
    rounded_delay = timestep * np.round(delay/timestep)
    print('Delay was rounded to ', rounded_delay)
    n_pad = int(rounded_delay/timestep)
    output_delta_p_profile = np.pad(delta_p_profile[0:-n_pad],[n_pad,0],constant_values=0)
    return output_delta_p_profile

## testing
if __name__ == "__main__":
    delta_p_profile = np.ones([100])
    delay = 4 
    timestep = 0.2

    print(simple_time_delay(delta_p_profile, delay, timestep))

