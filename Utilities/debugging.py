# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 08:43:51 2024

@author: felix
"""



import cantera as ct
import numpy as np
from CoolProp.CoolProp import PropsSI as psi

def set_combustion():
    comb = ct.Solution('gri30_WARR.yaml')
    return comb

def obtain_combustion_properties(comb, p):
    def heat_capacity_ratio(comb):
        return comb.cp_mass / comb.cv_mass
    fuel = 'C2H5OH'
    
    ROF = 1.4
    comb.Y = {fuel: 1, 'O2': ROF}
    comb.TP = 273, p
    comb.equilibrate('HP')
    
    delta_h_ox = get_delta_h(injec_t=100, ref_t=273, p=20e5, fluid='O2')
    delta_h_fuel = get_delta_h(injec_t=293, ref_t=273, p=20e5, fluid='C2H6O')
    
    deltah = ROF/(1+ROF)*delta_h_ox + 1/(1+ROF)*delta_h_fuel
    
    comb.HP = comb.h-deltah, p
    
    return comb.T, comb.mean_molecular_weight, heat_capacity_ratio(comb)

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

def calc_isentropic_exhaust_velocity(p_cc, T, M, gamma, p_ex=1e5):
    R_ideal = 8314
    # comb = set_combustion()
    v_ex = np.sqrt((2*gamma*R_ideal)/(gamma-1)*(T)/(M)*(1-(p_ex/p_cc)**((gamma-1)/gamma)))
    return v_ex

comb = set_combustion()

p = 20e5

temp, molar, gamma = obtain_combustion_properties(comb, p)

v_ex = calc_isentropic_exhaust_velocity(p, temp, molar, gamma, p_ex=1e5)

''' bei max chamber druck ist austrittsdruck als 0.4bar!!! '''

''' selbe combustion properties wenn ich bei CEA LOx und Ethanol als NICHT flüssig eingebe '''
''' wenn ich pc/p_ex = 20 wähle bei CEA und p_ex = 1 bar bei calc_isen... dann erhalte ich dieselbe geschwindigkeit '''