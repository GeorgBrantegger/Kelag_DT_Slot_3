import numpy as np
#based on https://en.wikipedia.org/wiki/PID_controller#Discrete_implementation


class P_controller_class:
    def __init__(self,setpoint,proportionality_constant):
        self.SP = setpoint
        self.Kp = proportionality_constant
        self.error_history = [] 
        self.control_variable = 0.1
        self.lower_limit = -0.1   # default
        self.upper_limit = +0.1   # default

    def set_control_variable_limits(self,lower_limit,upper_limit):
        self.lower_limit = lower_limit
        self.upper_limit = upper_limit

    def calculate_error(self,process_variable):
        self.error = self.SP-process_variable
        self.error_history.append(self.error)

    def get_control_variable(self):
        new_control = self.control_variable+self.Kp*(self.error_history[-1]-self.error_history[-2])
        if new_control < self.lower_limit:
            new_control = self.lower_limit

        if new_control > self.upper_limit:
            new_control = self.upper_limit

        self.control_variable = new_control
        # print(new_control)
        return new_control



