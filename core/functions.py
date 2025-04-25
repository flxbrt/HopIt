# -*- coding: utf-8 -*-
"""
Title: functions.py
Project: HopIt - Hopper Iterative System Analysis Tool
Author: @flxbrt
Version: 2.1

Description:
    Entry point for iterative hopper system analysis. Loads configuration,
    initializes combustion parameters, and runs a convergence loop to estimate
    total system mass through subsystem evaluations.
"""


import os
import numpy as np
import cantera as ct
import json
from CoolProp.CoolProp import PropsSI as psi


def load_system(config, direc, name="system"):
    path = os.path.join(direc, "config/")
    with open(path + name + '.json', 'r') as system_file:
        system = system_file.read()
    system = json.loads(system)
    system['Config'] = config
    return system

def write_system(system, direc, file_name="system"):
    with open(direc + file_name + '.json', 'w') as system_file:
        json.dump(system, system_file, indent=4)

def print_system(system):
    for subsystem in system:
        if subsystem!='Config':
            print(subsystem)
            for component in system[subsystem]:
                print('    ' + component)
                for spec in system[subsystem][component]:
                    print(f'        {spec}: ' + str(system[subsystem][component][spec]['Value']))
                    # print(f'        {spec}: {system[subsystem][component][spec]["Value"]:.2f}')
                    # print('        %s spec: %.2f' % (spec, system[subsystem][component][spec]["Value"]))

def get_material_properties(material):
    mat_lower = material.lower() # put everything in lower case --> less sensitive to typos
    if mat_lower == 'aluminium':
    # https://www.engineeringtoolbox.com/properties-aluminum-pipe-d_1340.html
        # return 2700, 275e6 # density, yield strength
    # https://www.mtm-inc.com/custom-aluminum-asme-code-stamped-pressure-vessels.html
    # https://www.fergusonperf.com/the-perforating-process/material-information/specialized-aluminum/6061-aluminium-alloy/
        return 2700, 241e6
    elif mat_lower == 'carbon fibre':
        rho = 1800
        ys = 1000e6
        ys_corrected = ys/1.4
        # correction according to type IV tanks found on the web
        return rho, ys_corrected
    elif mat_lower == 'steel':
    # # https://www.engineeringtoolbox.com/young-modulus-d_417.html
    #     rho = 7930 # https://www.kloecknermetals.com/blog/what-is-the-density-of-stainless-steel/#:~:text=The%20density%20of%20grade%20304,at%207750%20g%2Fm3.
    #     ys = 205e6 # https://masteel.co.uk/news/what-304-stainless-steel/
        # assuming 304 steel
        # https://asm.matweb.com/search/SpecificMaterial.asp?bassnum=mq316a
        rho = 8000
        ys = 290e6
        return rho, ys
    else:
        raise NameError(f'Material "{material}" not defined')

def set_combustion():
    comb = ct.Solution('gri30_WARR.yaml')
    return comb

def obtain_combustion_properties(system, comb, p):
    def heat_capacity_ratio(comb):
        return comb.cp_mass / comb.cv_mass
    
    def get_delta_h(injec_t, ref_t, p, fluid):
        
        vap_t = psi('T','P',p,'Q',1,fluid)
        # ref_t: temperature with which cantera object is initialized
        # vap_t: vaporization temperature of fluid at given pressure
        # injec_t: propellant temperature when injected / tank propellant temperature
        # assumption: propellants are always liquid --> vaporisation temperature is always higher than injection temperature
        # 3 cases need to be distinguished
        # 1: injec_t < vap_t < ref_t
        if ref_t > vap_t:
            h_low = psi('H', 'P', p,'T', injec_t, fluid)
            h_high = psi('H', 'P', p,'T', ref_t, fluid)
            deltah = h_high - h_low
            # print(1)
        # 2: injec_t < ref_t < vap_t
        if ref_t < vap_t and ref_t > injec_t:
            deltah = psi('H', 'P', p,'T', ref_t, fluid) - psi('H', 'P', p,'T', injec_t, fluid) + psi('H','P',p,'Q',1,fluid) - psi('H','P',p,'Q',0,fluid)
            # print(2)
        # 3: ref_t < injec_t < vap_t
        if ref_t < injec_t:
            # negative sign for the first paranthesis, since deltah is subtracted and therefore negative frist term yields less enthalpy being subtracted from cantera object
            deltah =  - (psi('H', 'P', p,'T', injec_t, fluid) - psi('H', 'P', p,'T', ref_t, fluid)) + psi('H','P',p,'Q',1,fluid) - psi('H','P',p,'Q',0,fluid)
            # print(3)
        return deltah 
    
    ref_t = 273
    
    fuel = system['Config']['Fuel']
    fuel_cp = system['Config']['Fuel Coolprop']
    ox = system['Config']['Oxidizer']
    ROF = system['Propulsion']['Thrust Chamber']['Oxidizer Fuel Ratio']['Value']
    
    injec_t_fuel = system['Propulsion']['Fuel']['Temperature']['Value']
    injec_t_ox = system['Propulsion']['Oxidizer']['Temperature']['Value']

    comb.Y = {fuel: 1, ox: ROF}
    comb.TP = ref_t, p
    comb.equilibrate('HP')
    
    delta_h_ox = get_delta_h(injec_t=injec_t_ox, ref_t=ref_t, p=p, fluid=ox)
    delta_h_fuel = get_delta_h(injec_t=injec_t_fuel, ref_t=ref_t, p=p, fluid=fuel_cp)
    
    deltah = ROF/(1+ROF)*delta_h_ox + 1/(1+ROF)*delta_h_fuel
    
    comb.HP = comb.h-deltah, p
    
    return comb.T, comb.mean_molecular_weight, heat_capacity_ratio(comb)

def get_system_mass(system, ignore=list()):
    mass = 0
    for subsystem in system.keys():
        if subsystem != 'Config':
            for component in system[subsystem].keys():
                if component in ignore:
                    continue
                mass += system[subsystem][component]['Mass']['Value']
    return mass

def calculate_ullage(system):
    fuel_coolprop = system['Config']['Fuel Coolprop']
    
    # self pressurized
    # TODO: does only work for self pressurizing pressurlant --> if p > p_amb at the temperature that is supposed to be in the tank
    rho_prop_init = psi('D','T',system['Propulsion']['Fuel']['Temperature']['Value'],'Q',0,fuel_coolprop) # Q denotes amount of vapor --> here: no vapor --> 1 would be 100% vapor
    # pressurized through nitrogen
    rho_prop_end = psi('D','P',system['Propulsion']['Fuel']['Pressure']['Value'],'T',system['Propulsion']['Fuel']['Temperature']['Value'],fuel_coolprop)

    ullage = system['Propulsion']['Fuel Tank']['Ullage']['Value'] # residual volume that is not occupied though propellant when tank is filled

    ullage_fuel = 1-(1-ullage)*rho_prop_init/rho_prop_end
    
    system['Propulsion']['Fuel Tank']['Ullage'] = {'Type': 'Output', 'Value': ullage_fuel}
    
    ox = system['Config']['Oxidizer']
    ox_press = system['Config']['Oxidizer Self Pressurised']
    
    if ox_press:
        p_ox = psi('P','T',system['Propulsion']['Oxidizer']['Temperature']['Value'],'Q',0,ox)
        system['Propulsion']['Oxidizer']['Pressure'] = {'Type': 'Output', 'Value': p_ox + 0.1e5} # 0.1bar extra required, otherwise coolprop fails
        
    else:
        # self pressurized
        rho_ox_init = psi('D','T',system['Propulsion']['Oxidizer']['Temperature']['Value'],'Q',0,ox) # Q denotes amount of vapor --> here: no vapor --> 1 would be 100% vapor
        # pressurized through nitrogen
        rho_ox_end = psi('D','P',system['Propulsion']['Oxidizer']['Pressure']['Value'],'T',system['Propulsion']['Oxidizer']['Temperature']['Value'],ox)
    
        ullage = system['Propulsion']['Oxidizer Tank']['Ullage']['Value'] # residual volume that is not occupied though propellant when tank is filled
    
        ullage_ox = 1-(1-ullage)*rho_ox_init/rho_ox_end
    
        system['Propulsion']['Oxidizer Tank']['Ullage'] = {'Type': 'Output', 'Value': ullage_ox}
    
    return system

# TODO: combustion efficiency richtig einbeziehen
def calc_isentropic_exhaust_velocity(system, comb, indicator, p_ex=1e5, eta=1):
    if indicator == 'max':
        p_cc = system['Propulsion']['Thrust Chamber']['Chamber Pressure']['Value']
    elif indicator == 'min':
        p_cc = system['Propulsion']['Thrust Chamber']['Chamber Pressure Min']['Value']
    eta = system['Propulsion']['Thrust Chamber']['Efficiency']['Value'] 
    R_ideal = 8314
    # comb = set_combustion()
    T, M, gamma = obtain_combustion_properties(system, comb, p_cc)
    v_ex = np.sqrt((2*gamma*R_ideal)/(gamma-1)*(T)/(M)*(1-(p_ex/p_cc)**((gamma-1)/gamma)))
    system['Propulsion']['Thrust Chamber'][f'Isp_{indicator}'] = {'Type': 'Output', 'Value': v_ex*eta}
    return system

def calc_effective_exhaust_velocity(system, comb, indicator):
    if indicator == 'max':
        p_cc = system['Propulsion']['Thrust Chamber']['Chamber Pressure']['Value']
        m_dot = system['Propulsion']['Thrust Chamber']['Mass Flow Max']['Value']
    elif indicator == 'min':
        p_cc = system['Propulsion']['Thrust Chamber']['Chamber Pressure Min']['Value']
        # !!! m_dot = system['Propulsion']['Thrust Chamber']['Mass Flow Max']['Value']
        m_dot = system['Propulsion']['Thrust Chamber']['Mass Flow Min']['Value']
    
    d_th = system['Propulsion']['Thrust Chamber']['Throat Diameter']['Value']
    epsilon = system['Propulsion']['Thrust Chamber']['Geometric Expansion Ratio']['Value']
    A_exit = d_th**2/4*np.pi*epsilon
    
    p_ratio = system['Propulsion']['Thrust Chamber']['Pressure Expansion Ratio']['Value']
    p_ex = p_ratio*p_cc
    
    R_ideal = 8314
    T, M, gamma = obtain_combustion_properties(system, comb, p_cc)
    v_ex = np.sqrt((2*gamma*R_ideal)/(gamma-1)*(T)/(M)*(1-(p_ex/p_cc)**((gamma-1)/gamma)))

    p_amb = 1e5
    
    v_eff = v_ex + A_exit/m_dot*(p_ex-p_amb)
    
    if "Efficiency" in system['Propulsion']['Thrust Chamber'].keys():
        efficiency = system['Propulsion']['Thrust Chamber']['Efficiency']['Value']
    else:
        efficiency = 1
        
    system['Propulsion']['Thrust Chamber'][f'Isp_{indicator}_eff'] = {'Type': 'Output', 'Value': v_eff*efficiency}
    return system

def ziolkowsky(system, mode='design'):
    if mode == 'design':
        if 'Isp_mean_eff' in system['Propulsion']['Thrust Chamber'].keys():
            v_ex = system['Propulsion']['Thrust Chamber']['Isp_mean_eff']['Value']
        else:
            v_ex = system['Propulsion']['Thrust Chamber']['Isp_max']['Value']
        
        flight_time = system['Config']['Flight Time [s]']
        # print(flight_time)
        deltav = 9.81 * flight_time
        m_dry = get_system_mass(system, ["Fuel", "Oxidizer"]) # pressurant does not leave the system!!!
        
        m_prop = m_dry*(-1+np.exp(deltav/v_ex)) # without m_pressurant
        
        ROF = system['Propulsion']['Thrust Chamber']['Oxidizer Fuel Ratio']['Value']

        m_fuel = 1 / (1+ROF) * m_prop
        m_ox = ROF / (1+ROF) * m_prop
        system['Propulsion']['Fuel']['Mass'] = {'Type': 'Output', 'Value': m_fuel}
        system['Propulsion']['Oxidizer']['Mass'] = {'Type': 'Output', 'Value': m_ox}
        
        return system
    
def solve_structure_frame(system):
    mat = system['Structure']['Frame']['Material']['Value']
    rho, ys = get_material_properties(mat)  
    safety = system['Structure']['Frame']['Safety Factor']['Value']
    
    F = system['Propulsion']['Thrust Chamber']['Thrust max calculated']['Value']
    
    # tank volumina
    vol_press = system['Propulsion']['Pressurant Tank']['Volume']['Value']
    vol_fuel = system['Propulsion']['Fuel Tank']['Volume']['Value']
    vol_ox = system['Propulsion']['Oxidizer Tank']['Volume']['Value']
    vol_structure = vol_press + vol_fuel + vol_ox
    if system['EPS']['Battery']['Mass'] != 0:
        vol_battery = system['EPS']['Battery']['Volume']['Value']
        vol_structure += vol_battery
    
    # account for ullage volume --> 1.2 rough approximation
    vol_structure *= 1.2
    r_inner = max(system['Propulsion']['Pressurant Tank']['Radius inner']['Value'],
                  system['Propulsion']['Oxidizer Tank']['Radius inner']['Value'],
                  system['Propulsion']['Fuel Tank']['Radius inner']['Value'])
    l_structure = vol_structure/(np.pi*r_inner**2)
    
    # has to be adapted according to TVC max gimabl angle --> currently assumed to be 30°
    Q = F*np.sin(30/180*np.pi)
    # assumed "feste Einspannung"
    Mb = Q*l_structure
    
    # init of while loop
    sigma = 1
    r_outer = r_inner + 0.05
    while ys/sigma > safety:
        I = np.pi/64*16*(r_outer**4-r_inner**4)
        a_max = r_outer
        W = I/a_max
        sigma = Mb/W
        r_outer -= 0.5e-3
    # ensure at least 2mm wall thickness --> even this is consdered unrealistic low
    if r_outer - r_inner < 2e-3:
        r_outer = r_inner + 2e-3
    
    mass_structure = rho*l_structure*np.pi*(r_outer**2-r_inner**2)
    
    system['Structure']['Frame']['Radius inner'] = {'Type': 'Output', 'Value': r_inner}
    system['Structure']['Frame']['Radius outer'] = {'Type': 'Output', 'Value': r_outer}
    system['Structure']['Frame']['Length'] = {'Type': 'Output', 'Value': l_structure}
    system['Structure']['Frame']['Volume'] = {'Type': 'Output', 'Value': vol_structure}
    system['Structure']['Frame']['Max Tension'] = {'Type': 'Output', 'Value': sigma}
    system['Structure']['Frame']['Mass'] = {'Type': 'Output', 'Value': mass_structure}
    
    return system

def solve_thrust_chamber(system, comb):
    
    def expansion_ratio(p_e, p_cc, gamma):
        nominator = ((gamma-1)/2)*(2/(gamma+1))**((gamma+1)/(gamma-1))
        denominator = (p_e/p_cc)**(2/gamma)*(1-(p_e/p_cc)**((gamma-1)/gamma))
        return np.sqrt(nominator/denominator)
    
    def calc_throat_diameter(m_dot, gamma, M, T, p_cc):
        # theta = np.sqrt(gamma*(2/(gamma+1))**((gamma+1)/(2*(gamma-1))))
        theta = np.sqrt(gamma)*((gamma+1)/2)**(-(gamma+1)/(2*(gamma-1)))
        R_ideal = 8314
        # d_th = np.sqrt(m_dot/p_cc/np.pi*4/np.sqrt(M/(R_ideal*T))*theta)
        d_th = np.sqrt(4/np.pi*m_dot/p_cc*np.sqrt(T*R_ideal/M)/theta)
        return d_th
    
    def calc_chocked_mass_flow(d_th, gamma, M, T, p_cc):
        # theta = np.sqrt(gamma*(2/(gamma+1))**((gamma+1)/(2*(gamma-1))))
        theta = np.sqrt(gamma)*((gamma+1)/2)**(-(gamma+1)/(2*(gamma-1)))
        R_ideal = 8314
        m_dot = np.pi*d_th**2/4*p_cc/np.sqrt(T*R_ideal/M)*theta
        return m_dot
    
    def calc_isen_velocity(T, M, gamma, p_cc, p_ex=1e5):
        # copy of the function calc_isentropic_exhaust_velocity
        # however this function gets passed the combustion properties, while the other one is calling it
        # saves time to not setup a not combustion instance every time
        R_ideal = 8314
        eta = system['Propulsion']['Thrust Chamber']['Efficiency']['Value'] 
        v_ex = np.sqrt((2*gamma*R_ideal)/(gamma-1)*(T)/(M)*(1-(p_ex/p_cc)**((gamma-1)/gamma)))
        return v_ex*eta
    
    wet_mass = get_system_mass(system)
    dry_mass = get_system_mass(system, ['Fuel', 'Oxidizer', 'Pressurant'])    
    # required thrust
    min_thrust = dry_mass*9.81*system['Propulsion']['Thrust Chamber']['TW min required']['Value']
    max_thrust = wet_mass*9.81*system['Propulsion']['Thrust Chamber']['TW max required']['Value']
    if not "Throttling Ratio" in system['Config'].keys():
        throttling_ratio = min_thrust/max_thrust
    else:
        throttling_ratio = system['Config']['Throttling Ratio']
    
    p_cc_max = system['Propulsion']['Thrust Chamber']['Chamber Pressure']['Value']
    # summerfiel criterion --> important for min thrust value
    p_exit_min = system['Propulsion']['Thrust Chamber']['p exit min']['Value']
    p_amb = 1e5
    
    
    # reiterating over geometry until calculated thrust corresponds to required thrust

    # initial guesses
    # estimation of chamber pressure at lowest load point
    p_cc_min = p_cc_max * throttling_ratio
    A_e = 1e-4
    # m_dot_geometry = 1
    # m_dot_thrust = 0
    T_min = min_thrust + 100
    # while abs(m_dot_geometry - m_dot_thrust) > 1e-2:
        
    # comb = set_combustion()
    while abs(T_min - min_thrust)/min_thrust > 0.01:
        
        # max thrust
        T, M, gamma = obtain_combustion_properties(system, comb, p_cc_max)
        p_exit_max = p_exit_min*p_cc_max/p_cc_min
        v_ex_max = calc_isen_velocity(T, M, gamma, p_cc_max, p_exit_max)
        m_dot_max = (max_thrust - A_e*(p_exit_max - p_amb)) / v_ex_max
        d_th = calc_throat_diameter(m_dot_max, gamma, M, T, p_cc_max)
        
        # min thrust
        T, M, gamma = obtain_combustion_properties(system, comb, p_cc_min)
        epsilon = expansion_ratio(p_exit_min, p_cc_min, gamma)
        A_e = d_th**2/4*np.pi*epsilon
        v_ex_min = calc_isen_velocity(T, M, gamma, p_cc_min, p_exit_min)
        m_dot_min = calc_chocked_mass_flow(d_th, gamma, M, T, p_cc_min)
        T_min = m_dot_min*v_ex_min + A_e*(p_exit_min - p_amb)
        
        # iteration 
        p_cc_min = p_cc_min + 0.1e5*np.sign(min_thrust-T_min)
    
    T_min = m_dot_min*v_ex_min + A_e*(p_exit_min - p_amb)
    T_max = m_dot_max*v_ex_max + A_e*(p_exit_max - p_amb)
    
    # according to expansion ratio equation for geometric expansion ratio, pressure ratio and area ratio remain constant
    # independent of absolute values --> ratios are preserved
    p_exit_max = p_exit_min*p_cc_max/p_cc_min    
    
    system['Propulsion']['Thrust Chamber']['Chamber Pressure Min'] = {'Type': 'Output', 'Value': p_cc_min}
    system = calc_isentropic_exhaust_velocity(system, comb, 'min', p_exit_min)
    system = calc_isentropic_exhaust_velocity(system, comb, 'max', p_exit_max)
    Isp_mean = np.mean([system['Propulsion']['Thrust Chamber']['Isp_max']['Value'], system['Propulsion']['Thrust Chamber']['Isp_min']['Value']])
    system['Propulsion']['Thrust Chamber']['Isp_mean'] = {'Type': 'Output', 'Value': Isp_mean}
    system['Propulsion']['Thrust Chamber']['TW min calculated'] = {'Type': 'Output', 'Value': T_min/dry_mass/9.81}
    system['Propulsion']['Thrust Chamber']['TW max calculated'] = {'Type': 'Output', 'Value': T_max/wet_mass/9.81}
    system['Propulsion']['Thrust Chamber']['Thrust min calculated'] = {'Type': 'Output', 'Value': T_min}
    system['Propulsion']['Thrust Chamber']['Thrust max calculated'] = {'Type': 'Output', 'Value': T_max}
    system['Propulsion']['Thrust Chamber']['Thrust min required'] = {'Type': 'Output', 'Value': min_thrust}
    system['Propulsion']['Thrust Chamber']['Thrust max required'] = {'Type': 'Output', 'Value': max_thrust}
    system['Propulsion']['Thrust Chamber']['Geometric Expansion Ratio'] = {'Type': 'Output', 'Value': epsilon}
    system['Propulsion']['Thrust Chamber']['Pressure Expansion Ratio'] = {'Type': 'Output', 'Value': p_exit_min/p_cc_min}
    system['Propulsion']['Thrust Chamber']['Throat Diameter'] = {'Type': 'Output', 'Value': d_th}
    system['Propulsion']['Thrust Chamber']['Throttling Ratio required'] = {'Type': 'Output', 'Value': throttling_ratio}
    system['Propulsion']['Thrust Chamber']['Throttling Ratio calculated'] = {'Type': 'Output', 'Value': T_min/T_max}
    system['Propulsion']['Thrust Chamber']['Mass Flow Max'] = {'Type': 'Output', 'Value': m_dot_max}
    system['Propulsion']['Thrust Chamber']['Mass Flow Min'] = {'Type': 'Output', 'Value': m_dot_min}   
    system = calc_effective_exhaust_velocity(system, comb, 'min')
    system = calc_effective_exhaust_velocity(system, comb, 'max')
    Isp_mean_eff = np.mean([system['Propulsion']['Thrust Chamber']['Isp_max_eff']['Value'], system['Propulsion']['Thrust Chamber']['Isp_min_eff']['Value']])
    system['Propulsion']['Thrust Chamber']['Isp_mean_eff'] = {'Type': 'Output', 'Value': Isp_mean_eff}
    return system

def add_mass_margin(system):
    system_mass = get_system_mass(system, ['Pressurant', 'Oxidizer', 'Fuel'])
    margin_fac = system['Margin']['Margin']['Factor']['Value']
    margin_mass = system_mass*margin_fac
    system['Margin']['Margin']['Mass']['Value'] = margin_mass
    return system

def compute_valves(system):
    
    # 1) pressure reducers
    # assuming that highest Kv is required at the end of burn at full thrust --> lowest deltap across pressure reducer
    # equations from https://www.schweizer-fn.de/stroemung/kvwert/kvwert.php#def
    def get_Kv(Qn, rho_n, T1, deltap, p1, p2):
        if deltap<p1_end/2: # subcritical flow --> Ma<1
            Kv = Qn/514*np.sqrt(rho_n*T1/deltap/p2)    
        else: # supercritical flow --> Ma>1
            Kv = Qn/257/p1*np.sqrt(rho_n*T1)
        return Kv
    
    p1_init = system['Propulsion']['Pressurant']['Pressure']['Value']
    p1_end = system['Propulsion']['Pressurant']['Pressure End']['Value']
    T_1_init = system['Propulsion']['Pressurant']['Temperature Init']['Value']
    
    rho_n = psi('D','T',273,'P',101300,'N2')
    gamma = 1.4
    T_1 = T_1_init*(p1_end/p1_init)**((gamma-1)/gamma) # end temperature given isentropic expansion in tank
    
    # fuel
    fuel = system['Config']['Fuel']
    fuel_coolprop = system['Config']['Fuel Coolprop']
    
    p2 = system['Propulsion']['Fuel']['Pressure']['Value']
    deltap = p1_end-p2
    # obtaining volume flow for nitrogen at standard conditions
    ROF = system['Propulsion']['Thrust Chamber']['Oxidizer Fuel Ratio']['Value']
    m_fuel = system['Propulsion']['Thrust Chamber']['Mass Flow Max']['Value']*1/(1+ROF)
    rho_fuel = psi('D','P',system['Propulsion']['Fuel']['Pressure']['Value'],'T',
                   system['Propulsion']['Fuel']['Temperature']['Value'],fuel_coolprop)   
    Q_fuel = m_fuel / rho_fuel # volume flow
    # Q_fuel must be equal to Q_nitrogen in propellant tank --> Q_nitrogen in pressurizer tank is higher
    # assuming isothermal state of change across pressure reducer
    rho_press = psi('D','T',T_1,'P',p2,'N2')
    Q_fuel_n = Q_fuel*rho_press/rho_n*3600 # m^3/h required
    
    Kv = get_Kv(Q_fuel_n, rho_n, T_1, deltap/1e5, p1_end/1e5, p2/1e5)
    system['Propulsion']['Pressure Reducer']['Kv Fuel'] = {'Type': 'Output', 'Value': Kv}
    
    # oxidizer
    ox = system['Config']['Oxidizer']
    p2 = system['Propulsion']['Oxidizer']['Pressure']['Value']
    deltap = p1_end-p2
    # obtaining volume flow for nitrogen at standard conditions
    m_ox = system['Propulsion']['Thrust Chamber']['Mass Flow Max']['Value']*ROF/(1+ROF)
    rho_ox = psi('D','P',system['Propulsion']['Oxidizer']['Pressure']['Value'],'T',
                   system['Propulsion']['Oxidizer']['Temperature']['Value'],ox)   
    Q_ox = m_ox / rho_ox # volume flow
    # Q_fuel must be equal to Q_nitrogen in propellant tank --> Q_nitrogen in pressurizer tank is higher
    # assuming isothermal state of change across pressure reducer
    Q_ox_n = Q_ox*rho_press/rho_n*3600 # m^3/h required
    
    Kv = get_Kv(Q_ox_n, rho_n, T_1, deltap/1e5, p1_end/1e5, p2/1e5)
    
    system['Propulsion']['Pressure Reducer']['Kv Oxidizer'] = {'Type': 'Output', 'Value': Kv}
    
    return system

def compute_pumps_and_batteries(system):
    
    def pump_power(prop):
        fuel = system['Config']['Fuel Coolprop']
        ox = system['Config']['Oxidizer']
        ROF = system['Propulsion']['Thrust Chamber']['Oxidizer Fuel Ratio']['Value']
        p_inlet = system['Propulsion'][prop]['Pressure']['Value']
        # assuming twice the chamber pressure as outlet pressure --> same assumption as for pressure fed system
        p_outlet = 2*system['Propulsion']['Thrust Chamber']['Chamber Pressure']['Value']
        # p_outlet = 80e5 + system['Propulsion']['Thrust Chamber']['Chamber Pressure']['Value']
        if prop == 'Fuel':
            m_dot = system['Propulsion']['Thrust Chamber']['Mass Flow Max']['Value']*1/(1+ROF)
            rho = psi('D','P',system['Propulsion']['Fuel']['Pressure']['Value'],'T',
                           system['Propulsion']['Fuel']['Temperature']['Value'],fuel)
        elif prop == 'Oxidizer':
            m_dot = system['Propulsion']['Thrust Chamber']['Mass Flow Max']['Value']*ROF/(1+ROF)
            rho = psi('D','P',system['Propulsion']['Oxidizer']['Pressure']['Value'],'T',
                           system['Propulsion']['Oxidizer']['Temperature']['Value'],ox)
        
        deltap = p_outlet-p_inlet
        
        return m_dot*deltap/rho
    
    def electrical_power(P, eta_pump, eta_motor=1):
        return P/eta_pump/eta_motor
    
    energy = 0
    P = 0
    mass_fuel_pump = 0
    mass_ox_pump = 0
    
    if system['Config']['Fuel Pump'] == True:
        eta_pump = system['Propulsion']['Fuel Pump']['Efficiency']['Value']
        P_fuel_pump = pump_power('Fuel')
        P_fuel_motor = electrical_power(P_fuel_pump, eta_pump)
        energy += P_fuel_motor*system['Config']['Flight Time [s]']# conservative approximation with full thrust --> pump pwer is approximated with full thrust
        P += P_fuel_motor
        # TODO: check if 500W fluid or motor power
        mass_fuel_pump = 2*np.sqrt(P_fuel_pump/500) # scaling law of sqrt() based on Lorenz Pump design --> approximate mass of 2kg @ 500W
        
        system['Propulsion']['Fuel Pump']['Mass'] = {'Type': 'Output', 'Value': mass_fuel_pump}
        system['Propulsion']['Fuel Pump']['Power Pump'] = {'Type': 'Output', 'Value': P_fuel_pump}
        system['Propulsion']['Fuel Pump']['Power Motor'] = {'Type': 'Output', 'Value': P_fuel_motor}
        
    if system['Config']['Oxidizer Pump'] == True:
        eta_pump = system['Propulsion']['Oxidizer Pump']['Efficiency']['Value']
        P_ox_pump = pump_power('Oxidizer')
        P_ox_motor = electrical_power(P_ox_pump, eta_pump)
        energy += P_ox_motor*system['Config']['Flight Time [s]']
        P += P_ox_motor
        # TODO: check if 500W fluid or motor power
        mass_ox_pump = 2*np.sqrt(P_ox_pump/500) # scaling law of sqrt() based on Lorenz Pump design --> approximate mass of 2kg @ 500W
        
        system['Propulsion']['Oxidizer Pump']['Mass'] = {'Type': 'Output', 'Value': mass_ox_pump}
        system['Propulsion']['Oxidizer Pump']['Power Pump'] = {'Type': 'Output', 'Value': P_ox_pump}
        system['Propulsion']['Oxidizer Pump']['Power Motor'] = {'Type': 'Output', 'Value': P_ox_motor}
    
    # worst case is considered
    margin = 0.25
    power_density = 1000 # W/kg
    energy_density = 180 # Wh/kg
    battery_mass_power = P/power_density*(1+margin)# mass approximation given max power
    energy = energy/3.6e3 # conversion into Wh 
    battery_mass_energy = energy/energy_density # mass approximation given max energy density
    
    battery_mass = max(battery_mass_energy, battery_mass_power)
    if battery_mass == 0:
        battery_mass = 2
    
    system['EPS']['Battery']['Mass'] = {'Type': 'Output', 'Value': battery_mass}
    
    return system