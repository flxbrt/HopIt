# -*- coding: utf-8 -*-
"""
Created on Fri Jan 26 10:55:28 2024

@author: felix
"""



#%%
from CoolProp.CoolProp import PropsSI as psi
import casadi as ca
from functions import get_material_properties

def solve_tanks(system):
    
    #%% calling required system parameters
    # fuel
    fuel_name = system['Config']['Fuel Coolprop']
    ox_name = system['Config']['Oxidizer']
    
    fuel_press = system['Config']['Fuel Self Pressurised']
    ox_press = system['Config']['Oxidizer Self Pressurised']
    if fuel_press:
        fuel_press = 0
    else:
        fuel_press = 1
    if ox_press:
        ox_press = 0
    else:
        ox_press = 1
        
    mat_fuel = system['Propulsion']['Fuel Tank']['Material']['Value']
    ullage_fuel = system['Propulsion']['Fuel Tank']['Ullage']['Value']
    ar_fuel = system['Propulsion']['Fuel Tank']['Aspect Ratio']['Value']
    mat_rho_fuel, ys_fuel = get_material_properties(mat_fuel)
    p_tank_fuel = system['Propulsion']['Fuel']['Pressure']['Value']
    sf_fuel = system['Propulsion']['Fuel Tank']['Safety Factor']['Value']
    m_fuel = system['Propulsion']['Fuel']['Mass']['Value']
    
    # oxidizer
    mat_ox = system['Propulsion']['Oxidizer Tank']['Material']['Value']
    ullage_ox = system['Propulsion']['Oxidizer Tank']['Ullage']['Value']
    ar_ox = system['Propulsion']['Oxidizer Tank']['Aspect Ratio']['Value']
    mat_rho_ox, ys_ox = get_material_properties(mat_ox)
    p_tank_ox = system['Propulsion']['Oxidizer']['Pressure']['Value']
    sf_ox = system['Propulsion']['Oxidizer Tank']['Safety Factor']['Value']
    m_ox = system['Propulsion']['Oxidizer']['Mass']['Value']
    
    # pressurizer - additionally, some precalculations
    mat_press = system['Propulsion']['Pressurant Tank']['Material']['Value']
    ar_press = system['Propulsion']['Pressurant Tank']['Aspect Ratio']['Value']
    mat_rho_press, ys_press = get_material_properties(mat_press)
    
    # TODO: doubles with p_press_init --> remove one
    p_tank_press = system['Propulsion']['Pressurant']['Pressure']['Value'] # initial pressure
    sf_press = system['Propulsion']['Pressurant Tank']['Safety Factor']['Value']

    p_press_init = system['Propulsion']['Pressurant']['Pressure']['Value'] # intial pressure
    p_press_end = system['Propulsion']['Pressurant']['Pressure End']['Value']  # end pressure
    T_press_init = system['Propulsion']['Pressurant']['Temperature Init']['Value'] # initial temperature
    # assuming that initial temperature in propellant tanks equivalent to initial temperature in pressurant tank
    # !not the case for LOx propellant tank! --> becomes less important for smaller ullage volume
    T_prop_tank_init = T_press_init 
    
    # TODO: doubles with p_tank_ox and p_tank_fuel --> remove one
    p_fuel = system['Propulsion']['Fuel']['Pressure']['Value']
    p_ox = system['Propulsion']['Oxidizer']['Pressure']['Value'] 
    # assuming isentropic state of change in pressurizer tank
    gamma = 1.4
    T_press_end = T_press_init*(p_press_end/p_press_init)**((gamma-1)/gamma)
    # assuming isenthalpic pressure reducer and ideal gas behavior across reducer
    T_prop_tank_end = T_press_end
    # calling coolprop for thermophysical properties of pressurizer
    # initial and final density in pressurizer tank 
    rho_press_init = psi('D','P',p_press_init,'T',T_press_init,'N2')
    rho_press_end = psi('D','P',p_press_end,'T',T_press_end,'N2')
    # initial and final density of pressurizer in propellant tank
    rho_fuel_prop_init = psi('D','P',p_fuel,'T',T_prop_tank_init,'N2')
    rho_fuel_prop_end = psi('D','P',p_fuel,'T',T_prop_tank_end,'N2')
    rho_ox_prop_init = psi('D','P',p_ox,'T',T_prop_tank_init,'N2')
    rho_ox_prop_end = psi('D','P',p_ox,'T',T_prop_tank_end,'N2')
    
    # to be called before opti is set up
    rho_fuel = psi('D','P',system['Propulsion']['Fuel']['Pressure']['Value'],'T',
                   system['Propulsion']['Fuel']['Temperature']['Value'],fuel_name)
    rho_ox = psi('D','P',system['Propulsion']['Oxidizer']['Pressure']['Value'],'T',
                   system['Propulsion']['Oxidizer']['Temperature']['Value'],ox_name)
    #%% setting up the root finder thorugh casadi opti stack
    
    opti = ca.Opti()
    s = opti.variable(30)
    
    # propane variables
    tank_vol_fuel = s[0]                    # fuel volume
    tank_vol_ullage_fuel = s[1]             # tank volume given fuel volume + ullage
    r_i_fuel = s[2]                         # inner radius of tank given tank volume and aspect ratio
    h_fuel = s[3]                           # height of cylindrical part of tank
    s_fuel = s[4]                           # thickness of tank based on Barlow and cylindrical part --> same thickness for sperical part
    r_o_fuel = s[5]                         # outer radius
    m_endcaps_fuel = s[6]                   # mass of endcaps / sperical parts of tank
    m_cyl_fuel = s[7]                       # mass of cylindrical part of tank
    m_tank_fuel = s[8]                      # entire structural tank mass
    
    # oxidizer variables
    tank_vol_ox = s[9]
    tank_vol_ullage_ox = s[10]
    r_i_ox = s[11]
    h_ox = s[12]
    s_ox = s[13]
    r_o_ox = s[14]
    m_endcaps_ox = s[15]
    m_cyl_ox = s[16]
    m_tank_ox = s[17]
    
    # pressurizer variables
    prop_mass_ullage_press_init = s[18]     # initial mass of pressurizer in the propellant tanks  --> ullage volume is filled
    prop_mass_press_end = s[19]             # final mass of pressurizer in the propellant tanks --> entirely filled
    delta_mass_press = s[20]                # pressurizer mass that goes form pressurizer tank to propellant tanks
    vol_press_tank = s[21]                  # volume of pressurizer tank
    mass_press = s[22]                      # mass of pressurizer
    
    r_i_press = s[23]                       # inner radius pressurizer tank
    h_press = s[24]                         # height of cylindrical part of tank
    s_press = s[25]                         # thickness of tank based on Barlow and cylindrical part --> same thickness for sperical part
    r_o_press = s[26]                       # outer radius
    m_endcaps_press = s[27]                 # mass of endcaps / sperical parts of tank
    m_cyl_press = s[28]                     # mass of cylindrical part of tank
    m_tank_press = s[29]                    # entire structural tank mass
    
    # propane equations
    opti.subject_to(tank_vol_fuel - m_fuel/rho_fuel == 0)
    opti.subject_to(tank_vol_ullage_fuel - tank_vol_fuel * (1+ullage_fuel) == 0)
    opti.subject_to(r_i_fuel - (tank_vol_ullage_fuel/(ca.pi*ar_fuel + 4/3*ca.pi))**(1/3) == 0)
    opti.subject_to(h_fuel - r_i_fuel*ar_fuel == 0)
    opti.subject_to(s_fuel - r_i_fuel / ys_fuel * p_tank_fuel * sf_fuel == 0)
    opti.subject_to(r_o_fuel - (s_fuel + r_i_fuel) == 0)
    opti.subject_to(m_endcaps_fuel - mat_rho_fuel * ca.pi*4/3 * (r_o_fuel**3 - r_i_fuel**3) == 0)
    opti.subject_to(m_cyl_fuel - mat_rho_fuel * h_fuel * ca.pi * (r_o_fuel**2 - r_i_fuel**2) == 0)
    opti.subject_to(m_tank_fuel - (m_endcaps_fuel + m_cyl_fuel) == 0)
    
    # oxidizer equations
    opti.subject_to(tank_vol_ox - m_ox/rho_ox == 0)
    opti.subject_to(tank_vol_ullage_ox -  tank_vol_ox * (1+ullage_ox) == 0)
    opti.subject_to(r_i_ox - (tank_vol_ullage_ox/(ca.pi*ar_ox + 4/3*ca.pi))**(1/3) == 0)
    opti.subject_to(h_ox - r_i_ox*ar_ox == 0)
    opti.subject_to(s_ox - r_i_ox / ys_ox * p_tank_ox * sf_ox == 0)
    opti.subject_to(r_o_ox - (s_ox + r_i_ox) == 0)
    opti.subject_to(m_endcaps_ox - mat_rho_ox * ca.pi*4/3 * (r_o_ox**3 - r_i_ox**3) == 0)
    opti.subject_to(m_cyl_ox - mat_rho_ox * h_ox * ca.pi * (r_o_ox**2 - r_i_ox**2) == 0)
    opti.subject_to(m_tank_ox - (m_endcaps_ox + m_cyl_ox) == 0)
    
    # pressurizer equations
    opti.subject_to(prop_mass_ullage_press_init - rho_fuel_prop_init*tank_vol_fuel*ullage_fuel*fuel_press - rho_ox_prop_init*tank_vol_ox*ullage_ox*ox_press == 0)
    opti.subject_to(prop_mass_press_end - rho_fuel_prop_end*tank_vol_ullage_fuel*fuel_press - rho_ox_prop_end*tank_vol_ullage_ox*ox_press == 0)
    opti.subject_to(delta_mass_press - (prop_mass_press_end - prop_mass_ullage_press_init) == 0)
    opti.subject_to(vol_press_tank - delta_mass_press/(rho_press_init - rho_press_end) == 0)
    opti.subject_to(mass_press - vol_press_tank*rho_press_init - prop_mass_ullage_press_init == 0)
    
    opti.subject_to(r_i_press - (vol_press_tank/(ca.pi*ar_press + 4/3*ca.pi))**(1/3) == 0)
    opti.subject_to(h_press - r_i_press*ar_press == 0)
    opti.subject_to(s_press - r_i_press / ys_press * p_tank_press * sf_press == 0)
    opti.subject_to(r_o_press - (s_press + r_i_press) == 0)
    opti.subject_to(m_endcaps_press - mat_rho_press * ca.pi*4/3 * (r_o_press**3 - r_i_press**3) == 0)
    opti.subject_to(m_cyl_press - mat_rho_press * h_press * ca.pi * (r_o_press**2 - r_i_press**2) == 0)
    opti.subject_to(m_tank_press - (m_endcaps_press + m_cyl_press) == 0)
    
    # initialized for 120s --> works also for 10s and 1000s --> seems to be quite stable, if initialized
    # without initializatin it fails on the variable r_i_press --> I don't understand why :/
    # could instead also be initilized with system stored in .json file from a previous run
    # opti.set_initial(s[:25], [2.16882538e-02, 2.38570792e-02, 1.06237624e-01, 5.31188122e-01,
    #        8.49900996e-04, 1.07087525e-01, 2.18714277e-01, 5.44604381e-01,
    #        7.63318658e-01, 2.42391944e-02, 2.66631138e-02, 1.10249395e-01,
    #        5.51246976e-01, 3.20725513e-03, 1.13456650e-01, 1.36154608e+00,
    #        3.35483439e+00, 4.71638047e+00, 4.59274482e-03, 2.29872920e-01,
    #        3.71588038e+00, 3.48600746e+00, 2.73907534e-02, 9.20663982e+00, 0.11124331436378877])
    opti.set_initial(s[:24], [2.16882538e-02, 2.38570792e-02, 1.06237624e-01, 5.31188122e-01,
           8.49900996e-04, 1.07087525e-01, 2.18714277e-01, 5.44604381e-01,
           7.63318658e-01, 2.42391944e-02, 2.66631138e-02, 1.10249395e-01,
           5.51246976e-01, 3.20725513e-03, 1.13456650e-01, 1.36154608e+00,
           3.35483439e+00, 4.71638047e+00, 2.29872920e-01,
           3.71588038e+00, 3.48600746e+00, 2.73907534e-02, 9.20663982e+00, 0.11124331436378877])
    
    '''
    über simple if abfrage könnte man manche komponenzen constrainen und dann die entsprechenden gleichungen nicht übergeben
    init werte können sogar bleiben
    flugzeit darf dann aber kein constrint mehr sein, sondern muss ergebnis sein
    '''
    
    #%% solving the root finding problem
    opts = {'ipopt.print_level':0, 'print_time':0}
    opti.solver('ipopt', opts)
    sol = opti.solve()
    params = opti.value(s)   
    
    #%% writing the solution to the system
    
    system['Propulsion']['Fuel']['Volume'] = {'Type': 'Output', 'Value': params[0]}
    system['Propulsion']['Fuel Tank']['Volume'] = {'Type': 'Output', 'Value': params[1]}
    system['Propulsion']['Fuel Tank']['Radius inner'] = {'Type': 'Output', 'Value': params[2]}
    system['Propulsion']['Fuel Tank']['Height cylindrcal'] = {'Type': 'Output', 'Value': params[3]}
    system['Propulsion']['Fuel Tank']['Thickness'] = {'Type': 'Output', 'Value': params[4]}
    system['Propulsion']['Fuel Tank']['Radius outer'] = {'Type': 'Output', 'Value': params[5]}
    system['Propulsion']['Fuel Tank']['Mass Endcaps'] = {'Type': 'Output', 'Value': params[6]}
    system['Propulsion']['Fuel Tank']['Mass Cylindrical'] = {'Type': 'Output', 'Value': params[7]}
    system['Propulsion']['Fuel Tank']['Mass'] = {'Type': 'Output', 'Value': params[8]}
    
    system['Propulsion']['Oxidizer']['Volume'] = {'Type': 'Output', 'Value': params[9]}
    system['Propulsion']['Oxidizer Tank']['Volume'] = {'Type': 'Output', 'Value': params[10]}
    system['Propulsion']['Oxidizer Tank']['Radius inner'] = {'Type': 'Output', 'Value': params[11]}
    system['Propulsion']['Oxidizer Tank']['Height cylindrcal'] = {'Type': 'Output', 'Value': params[12]}
    system['Propulsion']['Oxidizer Tank']['Thickness'] = {'Type': 'Output', 'Value': params[13]}
    system['Propulsion']['Oxidizer Tank']['Radius outer'] = {'Type': 'Output', 'Value': params[14]}
    system['Propulsion']['Oxidizer Tank']['Mass Endcaps'] = {'Type': 'Output', 'Value': params[15]}
    system['Propulsion']['Oxidizer Tank']['Mass Cylindrical'] = {'Type': 'Output', 'Value': params[16]}
    system['Propulsion']['Oxidizer Tank']['Mass'] = {'Type': 'Output', 'Value': params[17]}    
    
    system['Propulsion']['Pressurant']['Mass Propellant Tank Init'] = {'Type': 'Output', 'Value': params[18]}
    system['Propulsion']['Pressurant']['Mass Propellant Tank End'] = {'Type': 'Output', 'Value': params[19]}
    system['Propulsion']['Pressurant']['Mass Change'] = {'Type': 'Output', 'Value': params[20]}
    system['Propulsion']['Pressurant']['Volume'] = {'Type': 'Output', 'Value': params[21]}
    system['Propulsion']['Pressurant Tank']['Volume'] = {'Type': 'Output', 'Value': params[21]}
    system['Propulsion']['Pressurant']['Mass'] = {'Type': 'Output', 'Value': params[22]}
    system['Propulsion']['Pressurant Tank']['Radius inner'] = {'Type': 'Output', 'Value': params[23]}
    system['Propulsion']['Pressurant Tank']['Height cylindrcal'] = {'Type': 'Output', 'Value': params[24]}
    system['Propulsion']['Pressurant Tank']['Thickness'] = {'Type': 'Output', 'Value': params[25]}
    system['Propulsion']['Pressurant Tank']['Radius outer'] = {'Type': 'Output', 'Value': params[26]}
    system['Propulsion']['Pressurant Tank']['Mass Endcaps'] = {'Type': 'Output', 'Value': params[27]}
    system['Propulsion']['Pressurant Tank']['Mass Cylindrical'] = {'Type': 'Output', 'Value': params[28]}
    system['Propulsion']['Pressurant Tank']['Mass'] = {'Type': 'Output', 'Value': params[29]}
    
    for medium in ['Fuel', 'Oxidizer', 'Pressurant']:
        if system['Propulsion'][f'{medium} Tank']['Thickness']['Value']/system['Propulsion'][f'{medium} Tank']['Radius inner']['Value']>0.05:
            print(f'Caution: Barlow for {medium} Tank not fulfilled --> Thickness / Radius more than 5%')
            print('%.2f percent' % (100*system['Propulsion'][f'{medium} Tank']['Thickness']['Value']/system['Propulsion'][f'{medium} Tank']['Radius inner']['Value']))
    
    # system['Propulsion']['Pressurant Tank']['Densities'] = {'rho_press_init': rho_press_init,
    #                                                         'rho_press_end': rho_press_end,
    #                                                         'rho_prop_init': rho_prop_init,
    #                                                         'rho_fuel_prop_init': rho_fuel_prop_init,
    #                                                         'rho_fuel_prop_end': rho_fuel_prop_end,
    #                                                         'rho_ox_prop_init': rho_ox_prop_init,
    #                                                         'rho_ox_prop_end': rho_ox_prop_end,
    #                                                         'rho_fuel': rho_fuel,
    #                                                         'rho_ox': rho_ox,
    #                                                         'T_press_end': T_press_end}    
    
    #%% exception for mass calculation of carbon fobre tank, since 
    # the mass prediction based on Barlow is having big uncertainties
    
    tanks = ['Fuel', 'Oxidizer', 'Pressurant']
    for tank in tanks:
        if system['Propulsion'][f'{tank} Tank']['Material']['Value'] == 'Carbon Fibre':
            performance_fac = system['Propulsion'][f'{tank} Tank']['Performance']['Value']
            m_tank = system['Propulsion'][f'{tank}']['Pressure']['Value']/1e5*system['Propulsion'][f'{tank} Tank']['Volume']['Value']*1000/performance_fac
            system['Propulsion'][f'{tank} Tank']['Mass']['Value'] = m_tank
        else:
            performance_fac = system['Propulsion'][f'{tank}']['Pressure']['Value']/1e5*system['Propulsion'][f'{tank} Tank']['Volume']['Value']*1000/system['Propulsion'][f'{tank} Tank']['Mass']['Value']
            system['Propulsion'][f'{tank} Tank']['Performance'] = {'Type': 'Output', 'Value': performance_fac}
    return system