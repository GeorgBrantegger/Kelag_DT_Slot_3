import numpy as np
#based on https://en.wikipedia.org/wiki/PID_controller#Discrete_implementation

def trap_int(vec,timestep):
    l = np.size(vec)
    int = 0
    for i in range(l-1):
        int = int + (vec[i]+vec[i+1])/2*timestep
    return int


def ISE_fun(error_history,timestep):
    # calcuate the integral of square error
    e = np.array(error_history)
    dt = timestep
    ise = trap_int(e**2,dt)
    return ise

def IAE_fun(error_history,timestep):
    # calcuate the integral of absolute error
    e = np.array(error_history)
    dt = timestep
    iae = trap_int(np.abs(e),dt)
    return iae    

def ITSE_fun(error_history,timestep):
    # calcuate the integral of time multiply square error
    e = np.array(error_history)
    dt = timestep
    n = np.size(e)
    t = np.arange(0,n)*dt
    itse = trap_int(t*e**2,dt)
    return itse

def ITAE_fun(error_history,timestep):
    # calcuate the integral of time multiply absolute error
    e = np.array(error_history)
    dt = timestep
    n = np.size(e)
    t = np.arange(0,n)*dt    
    itae = trap_int(np.abs(e),dt)
    return itae    

class P_controller_class:
    # def __init__(self,setpoint,proportionality_constant):
    #     self.SP = setpoint
    #     self.Kp = proportionality_constant
    #     self.error_history = [] 
    #     self.control_variable = 0.1
    #     self.lower_limit = -0.1   # default
    #     self.upper_limit = +0.1   # default

    # def set_control_variable_limits(self,lower_limit,upper_limit):
    #     self.lower_limit = lower_limit
    #     self.upper_limit = upper_limit

    # def calculate_error(self,process_variable):
    #     self.error = self.SP-process_variable
    #     self.error_history.append(self.error)

    # def get_control_variable(self):
    #     new_control = self.control_variable+self.Kp*(self.error_history[-1]-self.error_history[-2])
    #     if new_control < self.lower_limit:
    #         new_control = self.lower_limit

    #     if new_control > self.upper_limit:
    #         new_control = self.upper_limit

    #     self.control_variable = new_control
    #     # print(new_control)
    #     return new_control
    def __init__(self):
        pass


class PI_controller_class:
    def __init__(self,setpoint,deadband,proportionality_constant,Ti, timestep):
        self.SP = setpoint
        self.db = deadband
        self.Kp = proportionality_constant
        self.Ti = Ti
        self.dt = timestep
        self.error_history = [0] 

        self.cv_lower_limit = 0   # default
        self.cv_upper_limit = +1   # default


    def set_control_variable_limits(self,lower_limit,upper_limit):
        self.cv_lower_limit = lower_limit
        self.cv_upper_limit = upper_limit

    def calculate_error(self,process_variable):
        self.error = process_variable-self.SP
        self.error_history.append(self.error)

    def get_control_variable(self,process_variable):

        self.calculate_error(process_variable)

        cv = self.control_variable
        Kp = self.Kp
        Ti = self.Ti
        dt = self.dt

        e0 = self.error_history[-1]
        e1 = self.error_history[-2] 
        if abs(self.error) > self.db:
            new_control = cv+Kp*(e0-e1)+dt/Ti*e0
        else:
            new_control = cv

        if new_control < self.cv_lower_limit:
            new_control = self.cv_lower_limit

        if new_control > self.cv_upper_limit:
            new_control = self.cv_upper_limit

        self.control_variable = new_control
        return self.control_variable

    def get_performance_indicators(self,ISE=True,IAE=True,ITSE=True,ITAE=True):
        ise     = np.nan
        iae     = np.nan
        itse    = np.nan
        itae    = np.nan

        # self.error_history[1:] because the first value of the error history is set to [0]
        #   to avoid special case handling in the calculation of the controll variable
        if ISE == True:
            ise = ISE_fun(self.error_history[1:],self.dt)
        if IAE == True:
            iae = IAE_fun(self.error_history[1:],self.dt)
        if ITSE == True:
            itse = ITSE_fun(self.error_history[1:],self.dt)
        if ITAE == True:
            itae = ITAE_fun(self.error_history[1:],self.dt)

        return ise,iae,itse,itae

    
