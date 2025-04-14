# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 08:24:40 2024

@author: felix
"""



import numpy as np
import cantera as ct
from CoolProp.CoolProp import PropsSI as psi


def set_combustion():
    comb = ct.Solution('gri30_WARR.yaml')
    return comb

def obtain_combustion_properties(fuel_cp, fuel_ct, ox, ROF, comb, injec_t_fuel, injec_t_ox, p):
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
        
    comb.Y = {fuel_ct: 1, ox: ROF}
    comb.TP = ref_t, p
    comb.equilibrate('HP')
    
    delta_h_ox = get_delta_h(injec_t=injec_t_ox, ref_t=ref_t, p=p, fluid=ox)
    delta_h_fuel = get_delta_h(injec_t=injec_t_fuel, ref_t=ref_t, p=p, fluid=fuel_cp)
    
    deltah = ROF/(1+ROF)*delta_h_ox + 1/(1+ROF)*delta_h_fuel
    
    comb.HP = comb.h-deltah, p
    
    return comb.T, comb.mean_molecular_weight, heat_capacity_ratio(comb)

def calc_isentropic_exhaust_velocity(p_cc, T, M, gamma, p_ex=1e5):
    R_ideal = 8314
    # comb = set_combustion()
    v_ex = np.sqrt((2*gamma*R_ideal)/(gamma-1)*(T)/(M)*(1-(p_ex/p_cc)**((gamma-1)/gamma)))
    return v_ex


fuel_cp = 'C2H6O'
fuel_ct = 'C2H5OH'
ox = 'N2O'
ROF = 5
injec_t_fuel = 293
injec_t_ox = 293
p_cc = 20e5

comb = set_combustion()

T, M, gamma = obtain_combustion_properties(fuel_cp, fuel_ct, ox, ROF, comb, injec_t_fuel, injec_t_ox, p_cc)

v_ex = calc_isentropic_exhaust_velocity(p_cc, T, M, gamma, p_ex = 1e5)