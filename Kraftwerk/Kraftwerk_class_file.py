#importing Druckrohrleitung
import sys
import os
current = os.path.dirname(os.path.realpath('Main_Programm.ipynb'))
parent = os.path.dirname(current)
sys.path.append(parent)
from functions.pressure_conversion import pressure_conversion
from Turbinen.Turbinen_class_file import Francis_Turbine

class Kraftwerk_class:
    def __init__(self):
        self.turbines = []

    def add_turbine(self,turbine):
        self.turbines.append(turbine)

    def print_info(self):
        for turbine in self.turbines:
            turbine.get_info(full=True)
