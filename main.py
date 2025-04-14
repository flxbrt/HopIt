# -*- coding: utf-8 -*-
"""
Created on Fri Jan 26 08:51:12 2024

@author: felix
"""



import os
import time
direc = 'C:/Users/felix/Desktop/script_directory/'
name = 'system'#_n2o_tank_material_study'
direc_save = 'C:/Users/felix/Desktop/script_directory/'

file_name = 'system_test_bench.json'

os.chdir(direc)

import functions as fc
from casadi_core import solve_tanks

config = {
    'Flight Time [s]': 80,
    # 'Number Press Tanks': None,
    # 'Number Press Reducers': None,
    'Oxidizer Pump': False,
    'Fuel Pump': False,
    'Oxidizer': 'O2',
    'Fuel': 'C2H5OH',
    'Fuel Coolprop': 'C2H6O',
    'Fuel Self Pressurised': False,
    'Oxidizer Self Pressurised': False
    }

sys_print = False
store_system = False
# design or performance
# design: calculates entire system given performance requirements (mainly flight time)
# performance: calculates system specs given constraints, such as fixed tank colume / mass
mode = 'design' # design gives back feasible design or it diverges

system = fc.load_system(config, direc, name)

comb = fc.set_combustion()

#%% main iteration loop

# calculate necessary ullage volume --> ullage factor
system = fc.calculate_ullage(system)

counter = 0

system_mass_prev = 1
system_mass = fc.get_system_mass(system)

# at first, exhaust velocity needs to be evaluated
system = fc.calc_isentropic_exhaust_velocity(system, comb, 'max')
threshold = 0.1 # [kg]

t_start = time.time()
while abs(system_mass_prev - system_mass) > threshold:
    system_mass_prev = system_mass
    # 1) ziolkowsy prop mass calculation    
    system = fc.ziolkowsky(system, mode)
    # 2) Tank mass calculation
    system = solve_tanks(system)
    # 3) Thrust Chamber calculation (mass, thrust, isp --> mean isp!!!)
    system = fc.solve_thrust_chamber(system, comb)
    # 4) Battery calculation
    # TODO: volume needs to be calculated
    system = fc.compute_pumps_and_batteries(system)
    # 5) structure calculation --> thrust chmaber needs to be solved before since strcture required max thrust
    system = fc.solve_structure_frame(system)
    # 6) add mass margin based on dry mass only
    system = fc.add_mass_margin(system)
    
    system_mass = fc.get_system_mass(system)
    counter += 1
    print(f'Iteration {counter} done')
    print(f'Relative deviation: {(system_mass-system_mass_prev)/system_mass_prev*100:.3f} [%]\n')
    
system = fc.compute_valves(system)

t_end = time.time()
   
print(f'Converged in {t_end-t_start:.2f} [s] and {counter} Iterations')


#%% print system

if sys_print:
    fc.print_system(system)
        
#%% store system

if store_system:
    fc.write_system(system, direc_save, file_name)

# temp
# print(fc.get_system_mass(system))
# print(system['Propulsion']['Oxidizer Tank']['Mass']['Value'])
# print(system['Propulsion']['Fuel Tank']['Mass']['Value'])
# print(system['Propulsion']['Pressurant Tank']['Mass']['Value'])
# print(system['Propulsion']['Oxidizer']['Mass']['Value'])
# print(system['Propulsion']['Fuel']['Mass']['Value'])
# print(system['Propulsion']['Pressurant']['Mass']['Value'])
# print(system['Propulsion']['Thrust Chamber']['Isp_mean_eff']['Value'])