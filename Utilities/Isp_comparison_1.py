# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 08:12:54 2024

@author: felix
"""



import cantera as ct
import numpy as np
import matplotlib.pyplot as plt



def heat_capacity_ratio(comb):
    return comb.cp_mass / comb.cv_mass

def calc_isentropic_exhaust_velocity(comb, p_cc, p_ex=1e5):
    R_ideal = 8314
    gamma = heat_capacity_ratio(comb)
    T = comb.T
    M = comb.mean_molecular_weight
    v_ex = np.sqrt((2*gamma*R_ideal)/(gamma-1)*(T)/(M)*(1-(p_ex/p_cc)**((gamma-1)/gamma)))
    return v_ex


fuel = 'C2H5OH'
ox_list = ['O2', 'N2O']
p_cc = 20e5 

ROF_list = np.linspace(1,6,100)
ex_vel = np.zeros((len(ox_list), len(ROF_list)))
    
for index_ox, ox in enumerate(ox_list):

    for index_ROF, ROF in enumerate(ROF_list):
        print(f'{ROF=}')
        comb = ct.Solution('gri30_WARR.yaml')
        comb.Y = {fuel: 1, ox: ROF}
        comb.TP = 293, p_cc
        comb.equilibrate('HP')
        ex_vel[index_ox, index_ROF] = calc_isentropic_exhaust_velocity(comb, p_cc)
    
    plt.plot(ROF_list, ex_vel[index_ox], label=ox)
    plt.xlabel('ROF')
    plt.ylabel('Exhaust Velocity')
    plt.grid()
    plt.legend()

#%%

index = 0
print('Max Isp of %.2f at ROF of %.2f' % (ex_vel[index, np.where(ex_vel[index,:]==max(ex_vel[index,:]))[0][0]], ROF_list[np.where(ex_vel[index,:]==max(ex_vel[index,:]))[0][0]]))


