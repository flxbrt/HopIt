# -*- coding: utf-8 -*-
"""
Created on Mon Jul  8 16:54:40 2024

@author: felix
"""



from CoolProp.CoolProp import PropsSI as psi

kappa = 1.4

vol_fuel = 0.04589931111271665
vol_ox = 0.03644340567211145

ullage_fuel = 0.10393815766725523
ullage_ox = 0.10627205541435047

vol_ullage_fuel = vol_fuel*ullage_fuel
vol_ullage_ox = vol_ox*ullage_ox

vol_fuel_tank = vol_fuel + vol_ullage_fuel
vol_ox_tank = vol_ox + vol_ullage_ox

T_init = 293
p_init = 300e5
p_end = 80e5

p_reducer_fuel = 40e5
p_reducer_ox = 30e5

T_end = T_init*(p_end/p_init)**((kappa-1)/kappa)

#%% using real gas data from coolprop

rho_press_init_fueltank = psi('D','P',p_reducer_fuel,'T',T_init,'N2')
rho_press_init_oxtank = psi('D','P',p_reducer_ox,'T',T_init,'N2')
# rho_press_init_proptank = p_reducer/(8314/28)/T_init

m_press_init_proptank = vol_ullage_fuel*rho_press_init_fueltank + vol_ullage_ox*rho_press_init_oxtank

rho_press_end_fueltank = psi('D','P',p_reducer_fuel,'T',T_end,'N2')
rho_press_end_oxtank = psi('D','P',p_reducer_ox,'T',T_end,'N2')

# rho_press_init_proptank = p_reducer/(8314/28)/T_end

m_press_end_proptank = vol_fuel_tank*rho_press_end_fueltank + vol_ox_tank*rho_press_end_oxtank

deltam = m_press_end_proptank - m_press_init_proptank

rho_presstank_init = psi('D','P',p_init,'T',T_init,'N2')
# rho_presstank_init = p_init/(8314/28)/T_init

rho_presstank_end = psi('D','P',p_end,'T',T_end,'N2')
# rho_presstank_end = p_end/(8314/28)/T_end

vol_press_senity = deltam/(rho_presstank_init - rho_presstank_end)

mass_press_senity = m_press_init_proptank + rho_presstank_init*vol_press_senity

#%% assuming ideal gas behavior all along
rho_press_init_fueltank_is = p_reducer_fuel/(8314/28)/T_init
rho_press_init_oxtank_is = p_reducer_ox/(8314/28)/T_init

m_press_init_proptank_is = vol_ullage_fuel*rho_press_init_fueltank_is + vol_ullage_ox*rho_press_init_oxtank_is

rho_press_end_fueltank_is = p_reducer_fuel/(8314/28)/T_end
rho_press_end_oxtank_is = p_reducer_ox/(8314/28)/T_end

m_press_end_proptank_is = vol_fuel_tank*rho_press_end_fueltank_is + vol_ox_tank*rho_press_end_oxtank_is

deltam_is = m_press_end_proptank_is - m_press_init_proptank_is

rho_presstank_init_is = p_init/(8314/28)/T_init

rho_presstank_end_is = p_end/(8314/28)/T_end

vol_press_senity_is = deltam_is/(rho_presstank_init_is - rho_presstank_end_is)

mass_press_senity_is = m_press_init_proptank_is + rho_presstank_init_is*vol_press_senity_is

#%%

print('Required Pressuriser Mass [kg]')
print(f'Real Gas: {mass_press_senity:.2f} | Ideal Gas: {mass_press_senity_is:.2f}')
print('Required Pressuriser Tanl Volume [l]')
print(f'Real Gas: {vol_press_senity*1000:.2f} | Ideal Gas: {vol_press_senity_is*1000:.2f}')
print('Density Pressurant Tank Init [kg/m^3]')
print(f'Real Gas: {rho_presstank_init:.2f} | Ideal Gas: {rho_presstank_init_is:.2f}')
print('Density Pressurant Tank End [kg/m^3]')
print(f'Real Gas: {rho_presstank_end:.2f} | Ideal Gas: {rho_presstank_end_is:.2f}')
