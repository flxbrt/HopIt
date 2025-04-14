# -*- coding: utf-8 -*-
"""
Created on Fri Aug  2 14:22:16 2024

@author: felix
"""



import numpy as np
import cantera as ct
from scipy.optimize import fsolve
from CoolProp.CoolProp import PropsSI as psi



def set_combustion():
    comb = ct.Solution('gri30_WARR.yaml')
    return comb

def obtain_combustion_properties(comb, p, fuel, fuel_cp, ROF, injec_t_fuel, injec_t_ox):
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
    

    comb.Y = {fuel: 1, 'O2': ROF}
    comb.TP = ref_t, p
    comb.equilibrate('HP')
    
    delta_h_ox = get_delta_h(injec_t=injec_t_ox, ref_t=ref_t, p=p, fluid='O2')
    delta_h_fuel = get_delta_h(injec_t=injec_t_fuel, ref_t=ref_t, p=p, fluid=fuel_cp)
    
    deltah = ROF/(1+ROF)*delta_h_ox + 1/(1+ROF)*delta_h_fuel
    
    comb.HP = comb.h-deltah, p
    
    return comb.T, comb.mean_molecular_weight, heat_capacity_ratio(comb)

def calc_effective_exhaust_velocity(p_cc, m_dot, d_th, epsilon, gamma, T, M):
    def expansion_ratio(p_e, p_cc, gamma):
        nominator = ((gamma-1)/2)*(2/(gamma+1))**((gamma+1)/(gamma-1))
        denominator = (p_e/p_cc)**(2/gamma)*(1-(p_e/p_cc)**((gamma-1)/gamma))
        return np.sqrt(nominator/denominator)
    
    def pressure_ratio(epsilon, gamma):
        def f(p_ratio, epsilon, gamma):
            nominator = ((gamma-1)/2)*(2/(gamma+1))**((gamma+1)/(gamma-1))
            denominator = (p_ratio)**(2/gamma)*(1-(p_ratio)**((gamma-1)/gamma))
            return epsilon - np.sqrt(nominator/denominator)
        p_ratio = fsolve(f, 0.01, args=(epsilon, gamma))
        return p_ratio 
        
    
    
    
    p_ratio = pressure_ratio(epsilon, gamma)
    
    p_exit = p_cc*p_ratio
    
    A_exit = d_th**2/4*np.pi*epsilon
    
    R_ideal = 8314
    
    v_ex = np.sqrt((2*gamma*R_ideal)/(gamma-1)*(T)/(M)*(1-(p_exit/p_cc)**((gamma-1)/gamma)))

    p_amb = 1e5
    
    v_eff = v_ex + A_exit/m_dot*(p_exit-p_amb)
    
    return v_eff, p_exit


#%%

m_dot_min = 0.627518891381022
m_dot_max = 1.3489163637096921
p_cc_min = 913340.0647136741
p_cc_max = 20e5
d_th = 0.037129838225282145

# annahme: lineare korrealtion zwischen massenstrom und druck --> gute näherung
epsilon = np.linspace(2, 4, 50)
m_dot = np.linspace(m_dot_min, m_dot_max, 10)
p_cc = np.linspace(p_cc_min, p_cc_max, 10)

v_eff = np.zeros((len(m_dot), len(epsilon)))
p_exit = np.zeros((len(m_dot), len(epsilon)))

for ii in range(len(m_dot)):
    comb = set_combustion()
    fuel = 'C2H5OH'
    fuel_cp = 'C2H6O'
    ROF = 1.1
    injec_t_fuel = 293
    injec_t_ox = 110
    T, M, gamma = obtain_combustion_properties(comb, p_cc[ii], fuel, fuel_cp, ROF, injec_t_fuel, injec_t_ox)
    for jj in range(len(epsilon)):
        v_eff[ii,jj], p_exit[ii,jj] = calc_effective_exhaust_velocity(p_cc[ii], m_dot[ii], d_th, epsilon[jj], gamma, T, M)
    
#%%
    
import matplotlib.pyplot as plt




fig, ax1 = plt.subplots(figsize=(8, 8))
ax2 = ax1.twinx()
for ii in range(len(m_dot)):

    ax1.plot(epsilon, v_eff[ii,:], label=f'{ii}')
    ax2.plot(epsilon, p_exit[ii,:]/1e5, label=f'{ii}')
    
plt.legend()
plt.grid()

#%%

veff_mean = np.mean(v_eff, axis=0)
index = np.where(max(veff_mean)==veff_mean)[0][0]
print(f'Best mean Isp at expansion ratio {epsilon[index]}')

plt.plot(p_cc, p_exit[:,index]/1e5)
