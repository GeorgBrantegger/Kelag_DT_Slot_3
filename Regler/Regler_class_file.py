import numpy as np

#based on https://en.wikipedia.org/wiki/PID_controller#Discrete_implementation

# performance parameters for controllers
def trap_int(vec,timestep):
    # numerical integration via the trapeziod rule to calculate the performance parameters
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

# P controller
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

# PI controller
class PI_controller_class:
# init
    def __init__(self,setpoint,deadband,proportionality_constant,Ti,timestep,lower_limit=0.,upper_limit=1.):
        self.SP = setpoint
        self.db = deadband
        self.Kp = proportionality_constant
        self.Ti = Ti                        # ~integration time
        self.dt = timestep
        # use a list to be able to append more easily - will get converted to np.array when needed
        self.error_history = [0]            

        self.cv_lower_limit = lower_limit   # limits for the controll variable
        self.cv_upper_limit = upper_limit   # limits for the controll variable

# setter
    def set_setpoint(self,setpoint):
        self.SP = setpoint

    def set_control_variable(self,control_variable, display_warning=True):
        if display_warning == True:
            print('WARNING! You are setting the control variable of the PI controller manually! \
                Consider using the .update_controll_variable() method instead.')
        self.control_variable = control_variable

# getter
    def get_current_control_variable(self):
        return self.control_variable

    def get_error_history(self):
        return self.error_history[1:]

    def get_performance_indicators(self,ISE=True,IAE=True,ITSE=True,ITAE=True):
        # calculate and return the performance indicators of the error history
        ise     = np.nan
        iae     = np.nan
        itse    = np.nan
        itae    = np.nan

        # self.error_history[1:] because the first value of the error history is set to [0]
        #   to avoid special case handling in the calculation of the control variable
        if ISE  == True:
            ise = ISE_fun(self.error_history[1:],self.dt)
        if IAE  == True:
            iae = IAE_fun(self.error_history[1:],self.dt)
        if ITSE == True:
            itse = ITSE_fun(self.error_history[1:],self.dt)
        if ITAE == True:
            itae = ITAE_fun(self.error_history[1:],self.dt)

        return ise,iae,itse,itae

    def get_info(self):        
        new_line = '\n'
        # :<10 pads the self.value to be 10 characters wide
        print_str = (f"Controller has the following attributes: {new_line}" 
            f"----------------------------- {new_line}"
            f"Type                      =   PI Controller {new_line}"
            f"Setpoint                  =   {self.SP:<10} {new_line}"
            f"Deadband                  =   {self.db:<10} {new_line}"
            f"Proportionality constant  =   {self.Kp:<10} {new_line}"
            f"Integration time          =   {self.Ti:<10} [s] {new_line}"
            f"Current control variable  =   {round(self.control_variable,3):<10}  {new_line}"
            f"Lower limit CV            =   {self.cv_lower_limit:<10}  {new_line}"
            f"Upper limit CV            =   {self.cv_upper_limit:<10}  {new_line}"
            f"Simulation timestep       =   {self.dt:<10} [s] {new_line}"
            f"----------------------------- {new_line}")

        print(print_str)

# methods
    def calculate_error(self,process_variable):
        # calculate the error and expand the err history
        self.error = process_variable-self.SP
        self.error_history.append(self.error)

    def update_control_variable(self,process_variable):
        # calculate the current control variable and make sure it does not exceed the limits
        self.calculate_error(process_variable)

        # initialize some variables
        cv = self.control_variable
        Kp = self.Kp
        Ti = self.Ti
        dt = self.dt

        e0 = self.error_history[-1]
        e1 = self.error_history[-2] 

        # test if the error exceeds the deadband range
            # only if that is the case, change control variable
        if abs(self.error) > self.db:
            new_control = cv+Kp*(e0-e1)+dt/Ti*e0
            # ensure that the controll variable stays within the predefined limits
            if new_control < self.cv_lower_limit:
                new_control = self.cv_lower_limit
            if new_control > self.cv_upper_limit:
                new_control = self.cv_upper_limit
        else:
            new_control = cv

        # set the control variable attribute
        self.set_control_variable(new_control,display_warning=False)


    
