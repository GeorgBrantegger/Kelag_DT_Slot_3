def turbine_flux(p,LA,p_exp,cubic_coeff,quadratic_coeff,linear_coeff,const_coeff):
    return (p*1e-5)**p_exp*(cubic_coeff*LA**3+quadratic_coeff*LA**2+linear_coeff*LA+const_coeff)


class Francis_turbine_class:
    def __init__(self):
        pass

    def set_turbine_flux_parameters(self,p_exp,cubic_coeff,quadratic_coeff,linear_coeff,const_coeff):
        # extracted from the Muschelkurve of the Turbine and used to calculate the turbine flux for a given pressure
        self.p_exp              = p_exp 
        self.cubic_coeff        = cubic_coeff
        self.quadratic_coeff    = quadratic_coeff
        self.linear_coeff       = linear_coeff
        self.const_coeff        = const_coeff         

    def get_turbine_flux(self,pressure,Leitapparatöffnung):
        self.flux = turbine_flux(pressure,Leitapparatöffnung,self.p_exp,self.cubic_coeff,self.quadratic_coeff,self.linear_coeff,self.const_coeff)
        return self.flux