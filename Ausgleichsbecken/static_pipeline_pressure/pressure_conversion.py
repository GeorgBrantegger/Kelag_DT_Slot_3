# convert to Pa
def bar_to_pa(p):
    return p*1e5

def mWS_to_pa(p):
    return p*9.80665*1e3

def torr_to_pa(p):
    return p*133.322

def atm_to_pa(p):
    return p*101.325*1e3

def psi_to_pa(p):
    return p*6894.8

# convert from Pa
def pa_to_bar(p):
    return p*1e-5

def pa_to_mWS(p):
    return p*1/(9.80665*1e3)

def pa_to_torr(p):
    return p/133.322

def pa_to_atm(p):
    return p*1/(101.325*1e3)

 # converstion function

def pa_to_psi(p):
    return p/6894.8

def pressure_conversion(pressure, input_unit = 'bar', target_unit = 'Pa'):
    p = pressure
    if input_unit.lower()   == 'bar':
        p_pa = bar_to_pa(p)
    elif input_unit.lower() == 'mws':
        p_pa = mWS_to_pa(p)
    elif input_unit.lower() == 'torr':
        p_pa = torr_to_pa(p)
    elif input_unit.lower() == 'atm':
        p_pa = atm_to_pa(p)
    elif input_unit.lower() == 'psi':
        p_pa = psi_to_pa(p)
    elif input_unit.lower() == 'pa':
        p_pa = p
    else:
        raise Exception('Given input unit not recognised. \n Known units are: Pa, bar, mWs, Torr, atm, psi')
    
    if target_unit.lower()      == 'bar':
        return pa_to_bar(p_pa), target_unit
    elif target_unit.lower()    == 'mws':
        return pa_to_mWS(p_pa), target_unit
    elif target_unit.lower()    == 'torr':
        return pa_to_torr(p_pa), target_unit
    elif target_unit.lower()    == 'atm':
        return pa_to_atm(p_pa), target_unit
    elif target_unit.lower()    =='psi':
        return pa_to_psi(p_pa), target_unit
    elif target_unit.lower()    == 'pa':
        return p_pa, target_unit
    else:
        raise Exception('Given target unit not recognised. \n Known units are: Pa, bar, mWs, Torr, atm, psi')

# testing_pressure_conversion
if __name__ == '__main__':
    p = 1

    unit_dict = ['Pa','Bar','Torr','Atm','MWS','psi']

    for input_unit in unit_dict:
        for target_unit in unit_dict:
            converted_p = pressure_conversion(p,input_unit,target_unit)
            print(input_unit,target_unit)
            print(converted_p)