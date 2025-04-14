# -*- coding: utf-8 -*-
"""
Created on Mon Jul  8 15:46:52 2024

@author: felix
"""



import numpy as np
from CoolProp.CoolProp import PropsSI as psi
import matplotlib.pyplot as plt

'''
Namenskonvention
# _ig           ideal gas assumption
# _td           calculated as it is calculated in the sys ana tool --> temperature after reducer is equivalent to the pressurizer temperature at the end
# _td_w_ae      as before, but the temperature after the reducer is decreasing slowly with the adiabatic expansion
# _jt           as before, and also considering joult thompson across the reducer

'''


#%% Temperaturunterschied anschauen zwischen isenthalpe ZÄ ideales gas vs. berücksichtigung joule thompson

T_init = 293
p_init = 300e5
p_end = 80e5

p_reducer = 40e5

fluid = 'N2'

kappa = 1.4

T_end = T_init*(p_end/p_init)**((kappa-1)/kappa)

rho_init = psi('D','P',p_init,'T',T_init,fluid)
rho_end = psi('D','P',p_end,'T',T_end,fluid)

h_init = psi('H','P',p_init,'T',T_init,fluid)
h_end = psi('H','P',p_end,'T',T_end,fluid)

T_pr_init = psi('T','P',p_reducer,'H',h_init,fluid)
T_pr_end = psi('T','P',p_reducer,'H',h_end,fluid)

print(f'{T_init=}')
print(f'{T_end=}')
print(f'{T_pr_init=}')
print(f'{T_pr_end=}')




#%% 3 berechnung der gasentleerung entsprechend wie die tankauslegung ist --> also über coolprop

vol_press = 0.037347040397043035

vol_fuel = 0.04589931111271665
vol_ox = 0.03644340567211145

ullage_fuel = 0.10393815766725523
ullage_ox = 0.10627205541435047

vol_ullage_fuel = vol_fuel*ullage_fuel
vol_ullage_ox = vol_ox*ullage_ox

vol_fuel_tank = vol_fuel + vol_ullage_fuel
vol_ox_tank = vol_ox + vol_ullage_ox

press = 'N2'

T_init = 293
p_init = 300e5
p_end = 80e5

p_reducer_fuel = 40e5
p_reducer_ox = 30e5

plot = False

rho_press_init = psi('D','P',p_init,'T',T_init,press)
mass_press = rho_press_init*vol_press




m_dot_prop = 1
ROF = 1.1
m_dot_fuel = 1/(1+ROF)*m_dot_prop
m_dot_ox = ROF/(1+ROF)*m_dot_prop



T_fuel = 293
T_ox = 100
fuel = 'C2H6O'
ox = 'O2'
rho_fuel = psi('D','P',p_reducer_fuel,'T',T_fuel,fuel)
rho_ox = psi('D','P',p_reducer_ox,'T',T_ox,ox)

kappa = 1.4

T_end = T_init*(p_end/p_init)**((kappa-1)/kappa)

# calculation
vol_dot_fuel = m_dot_fuel / rho_fuel
vol_dot_ox = m_dot_ox / rho_ox

N = 1000
deltat = vol_fuel / vol_dot_fuel # same result for deltat obtained when vol_ox/vol_dot_ox is used
dt = deltat / N

p_current = 300e5

mass_press_array_td = np.zeros(N)
T_reducer_array_td = np.zeros(N)
T_presstank_array_td = np.zeros(N)
# rho_reducer_array = np.zeros(N)
p_press_array_td = np.zeros(N)




vol_dot_press_fuel = vol_dot_fuel
vol_dot_press_ox = vol_dot_ox




for ii in range(N):
    T_current = T_init*(p_current/p_init)**((kappa-1)/kappa)
    # T_reducer = T_current
    #!!! diese annahme ist wesentlich
    T_reducer = T_end
    # fuel
    rho_reducer_fuel = psi('D','P',p_reducer_fuel,'T',T_reducer,press)
    # ox
    rho_reducer_ox = psi('D','P',p_reducer_ox,'T',T_reducer,press)
    
    m_dot_reducer = vol_dot_press_fuel*rho_reducer_fuel + vol_dot_press_ox*rho_reducer_ox
    
    mass_press -= m_dot_reducer*dt
    
    rho_press = mass_press/vol_press
    
    p_current = psi('P','D',rho_press,'T',T_current,press)
    
    mass_press_array_td[ii] = mass_press
    T_reducer_array_td[ii] = T_reducer
    T_presstank_array_td[ii] = T_current
    # rho_reducer_array[ii] = rho_reducer
    p_press_array_td[ii] = p_current
    
if plot:

    plt.figure(figsize=(12,8))
    
    plt.subplot(221)
    plt.plot(np.linspace(0, deltat, N), mass_press_array_td)
    plt.xlabel('Time [s]')
    plt.ylabel('Mass Pressurizer [kg]')
    plt.grid()
    
    plt.subplot(222)
    plt.plot(np.linspace(0, deltat, N), T_reducer_array_td)
    plt.xlabel('Time [s]')
    plt.ylabel('Temperature after reducer [K]')
    plt.grid()
    
    # plt.subplot(223)
    # plt.plot(np.linspace(0, deltat, N), rho_reducer_array)
    # plt.xlabel('Time [s]')
    # plt.ylabel('Density after reducer [kg/m^3]')
    # plt.grid()
    
    plt.subplot(224)
    plt.plot(np.linspace(0, deltat, N), p_press_array_td/1e5)
    plt.xlabel('Time [s]')
    plt.ylabel('Pressurizer pressure [bar]')
    plt.grid()

# so kommt man sehr nahe an tankauslegng ran

#%% 4 so wie bei 3) aber berücksichtung, dass temperatur im pressurizer tank stetig abnimmt entsprechend der adaibaten entspannung

vol_press = 0.037347040397043035

vol_fuel = 0.04589931111271665
vol_ox = 0.03644340567211145

ullage_fuel = 0.10393815766725523
ullage_ox = 0.10627205541435047

vol_ullage_fuel = vol_fuel*ullage_fuel
vol_ullage_ox = vol_ox*ullage_ox

vol_fuel_tank = vol_fuel + vol_ullage_fuel
vol_ox_tank = vol_ox + vol_ullage_ox

press = 'N2'

T_init = 293
p_init = 300e5
p_end = 80e5

p_reducer_fuel = 40e5
p_reducer_ox = 30e5

plot = False

rho_press_init = psi('D','P',p_init,'T',T_init,press)
mass_press = rho_press_init*vol_press




m_dot_prop = 1
ROF = 1.1
m_dot_fuel = 1/(1+ROF)*m_dot_prop
m_dot_ox = ROF/(1+ROF)*m_dot_prop



T_fuel = 293
T_ox = 100
fuel = 'C2H6O'
ox = 'O2'
rho_fuel = psi('D','P',p_reducer_fuel,'T',T_fuel,fuel)
rho_ox = psi('D','P',p_reducer_ox,'T',T_ox,ox)

kappa = 1.4

T_end = T_init*(p_end/p_init)**((kappa-1)/kappa)

# calculation
vol_dot_fuel = m_dot_fuel / rho_fuel
vol_dot_ox = m_dot_ox / rho_ox

N = 1000
deltat = vol_fuel / vol_dot_fuel # same result for deltat obtained when vol_ox/vol_dot_ox is used
dt = deltat / N

p_current = 300e5

mass_press_array_td_w_ae = np.zeros(N)
T_reducer_array_td_w_ae = np.zeros(N)
T_presstank_array_td_w_ae = np.zeros(N)

# rho_reducer_array = np.zeros(N)
p_press_array_td_w_ae = np.zeros(N)




vol_dot_press_fuel = vol_dot_fuel
vol_dot_press_ox = vol_dot_ox




for ii in range(N):
    T_current = T_init*(p_current/p_init)**((kappa-1)/kappa)
    T_reducer = T_current
    # fuel
    rho_reducer_fuel = psi('D','P',p_reducer_fuel,'T',T_reducer,press)
    # ox
    rho_reducer_ox = psi('D','P',p_reducer_ox,'T',T_reducer,press)
    
    m_dot_reducer = vol_dot_press_fuel*rho_reducer_fuel + vol_dot_press_ox*rho_reducer_ox
    
    mass_press -= m_dot_reducer*dt
    
    rho_press = mass_press/vol_press
    
    p_current = psi('P','D',rho_press,'T',T_current,press)
    
    mass_press_array_td_w_ae[ii] = mass_press
    T_reducer_array_td_w_ae[ii] = T_reducer
    T_presstank_array_td_w_ae[ii] = T_current

    # rho_reducer_array[ii] = rho_reducer
    p_press_array_td_w_ae[ii] = p_current
    
if plot:

    plt.figure(figsize=(12,8))
    
    plt.subplot(221)
    plt.plot(np.linspace(0, deltat, N), mass_press_array_td_w_ae)
    plt.xlabel('Time [s]')
    plt.ylabel('Mass Pressurizer [kg]')
    plt.grid()
    
    plt.subplot(222)
    plt.plot(np.linspace(0, deltat, N), T_reducer_array_td_w_ae)
    plt.xlabel('Time [s]')
    plt.ylabel('Temperature after reducer [K]')
    plt.grid()
    
    # plt.subplot(223)
    # plt.plot(np.linspace(0, deltat, N), rho_reducer_array)
    # plt.xlabel('Time [s]')
    # plt.ylabel('Density after reducer [kg/m^3]')
    # plt.grid()
    
    plt.subplot(224)
    plt.plot(np.linspace(0, deltat, N), p_press_array_td_w_ae/1e5)
    plt.xlabel('Time [s]')
    plt.ylabel('Pressurizer pressure [bar]')
    plt.grid()

# so hat man mehr mass an pressurizer zur verfügung als man braucht

#%% 5 berücksichtigung joule thompson

vol_press = 0.037347040397043035

vol_fuel = 0.04589931111271665
vol_ox = 0.03644340567211145

ullage_fuel = 0.10393815766725523
ullage_ox = 0.10627205541435047

vol_ullage_fuel = vol_fuel*ullage_fuel
vol_ullage_ox = vol_ox*ullage_ox

vol_fuel_tank = vol_fuel + vol_ullage_fuel
vol_ox_tank = vol_ox + vol_ullage_ox

press = 'N2'

T_init = 293
p_init = 300e5
p_end = 80e5

p_reducer_fuel = 40e5
p_reducer_ox = 30e5

plot = False

rho_press_init = psi('D','P',p_init,'T',T_init,press)
mass_press = rho_press_init*vol_press




m_dot_prop = 1
ROF = 1.1
m_dot_fuel = 1/(1+ROF)*m_dot_prop
m_dot_ox = ROF/(1+ROF)*m_dot_prop



T_fuel = 293
T_ox = 100
fuel = 'C2H6O'
ox = 'O2'
rho_fuel = psi('D','P',p_reducer_fuel,'T',T_fuel,fuel)
rho_ox = psi('D','P',p_reducer_ox,'T',T_ox,ox)

kappa = 1.4

T_end = T_init*(p_end/p_init)**((kappa-1)/kappa)

# calculation
vol_dot_fuel = m_dot_fuel / rho_fuel
vol_dot_ox = m_dot_ox / rho_ox

N = 1000
deltat = vol_fuel / vol_dot_fuel # same result for deltat obtained when vol_ox/vol_dot_ox is used
dt = deltat / N

p_current = 300e5

mass_press_array_jt = np.zeros(N)
T_reducer_fuel_array_jt = np.zeros(N)
T_reducer_ox_array_jt = np.zeros(N)
T_presstank_array_jt = np.zeros(N)

# rho_reducer_array = np.zeros(N)
p_press_array_jt = np.zeros(N)




vol_dot_press_fuel = vol_dot_fuel
vol_dot_press_ox = vol_dot_ox




for ii in range(N):
    T_current = T_init*(p_current/p_init)**((kappa-1)/kappa)
    h_current = psi('H','P',p_current,'T',T_current,press)
    # fuel
    T_reducer_fuel = psi('T','P',p_reducer_fuel,'H',h_current,press)
    rho_reducer_fuel = psi('D','P',p_reducer_fuel,'H',h_current,press)
    # ox
    T_reducer_ox = psi('T','P',p_reducer_ox,'H',h_current,press)
    rho_reducer_ox = psi('D','P',p_reducer_ox,'H',h_current,press)
    
    m_dot_reducer = vol_dot_press_fuel*rho_reducer_fuel + vol_dot_press_ox*rho_reducer_ox
    
    mass_press -= m_dot_reducer*dt
    
    rho_press = mass_press/vol_press
    
    p_current = psi('P','D',rho_press,'T',T_current,press)
    
    mass_press_array_jt[ii] = mass_press
    T_reducer_fuel_array_jt[ii] = T_reducer_fuel
    T_reducer_ox_array_jt[ii] = T_reducer_ox
    T_presstank_array_jt[ii] = T_current

    # rho_reducer_array[ii] = rho_reducer
    p_press_array_jt[ii] = p_current
    
if plot:

    plt.figure(figsize=(12,8))
    
    plt.subplot(221)
    plt.plot(np.linspace(0, deltat, N), mass_press_array_jt)
    plt.xlabel('Time [s]')
    plt.ylabel('Mass Pressurizer [kg]')
    plt.grid()
    
    plt.subplot(222)
    plt.plot(np.linspace(0, deltat, N), T_reducer_fuel_array_jt, label='Fuel Side')
    plt.plot(np.linspace(0, deltat, N), T_reducer_ox_array_jt, label='Ox Side')
    plt.xlabel('Time [s]')
    plt.ylabel('Temperature after reducer [K]')
    plt.grid()
    
    # plt.subplot(223)
    # plt.plot(np.linspace(0, deltat, N), rho_reducer_array)
    # plt.xlabel('Time [s]')
    # plt.ylabel('Density after reducer [kg/m^3]')
    # plt.grid()
    
    plt.subplot(224)
    plt.plot(np.linspace(0, deltat, N), p_press_array_jt/1e5)
    plt.xlabel('Time [s]')
    plt.ylabel('Pressurizer pressure [bar]')
    plt.grid()


#%% 6 vergleich von allen

# _ig           ideal gas assumption
# _td           calculated as it is calculated in the sys ana tool
# _td_w_ae      as before, but the temperature after the reducer is decreasing slowly with the adaibatic expansion
# _jt           as before, and also considering joult thompson across the reducer

plt.figure(figsize=(18,6))

plt.subplot(131)
# plt.plot(np.linspace(0, deltat, N), mass_press_array_ig, label='Ideal Gas')
plt.plot(np.linspace(0, deltat, N), mass_press_array_td, label='Sys Ana Design')
plt.plot(np.linspace(0, deltat, N), mass_press_array_td_w_ae, label='Sys Ana Design & Adiabatic Expansion')
plt.plot(np.linspace(0, deltat, N), mass_press_array_jt, label='Real Gas / Joule Thompson')
plt.xlabel('Time [s]')
plt.ylabel('Mass Pressurizer [kg]')
plt.legend()
plt.grid()

plt.subplot(132)
# plt.plot(np.linspace(0, deltat, N), T_reducer_array_ig, label='Ideal Gas')
plt.plot(np.linspace(0, deltat, N), T_reducer_array_td, label='After Reducer Sys Ana Design')
plt.plot(np.linspace(0, deltat, N), T_reducer_array_td_w_ae, label='After Reducer Sys Ana Design & Adiabatic Expansion')
plt.plot(np.linspace(0, deltat, N), T_reducer_fuel_array_jt, label='After Reducer Real Gas Fuel Side')
plt.plot(np.linspace(0, deltat, N), T_reducer_ox_array_jt, label='After Reducer Real Gas Ox Side')
plt.plot(np.linspace(0, deltat, N), T_presstank_array_td, linestyle = 'dashed', label='Press Tank Sys Ana Design')
plt.plot(np.linspace(0, deltat, N), T_presstank_array_td_w_ae, linestyle = 'dashed', label='Press Tank Sys Ana Design & Adiabatic Expansion')
plt.plot(np.linspace(0, deltat, N), T_presstank_array_jt, linestyle = 'dashed', label='Press Tank Real Gas')
plt.xlabel('Time [s]')
plt.ylabel('Temperature [K]')
plt.legend()
plt.grid()

# plt.subplot(223)
# plt.plot(np.linspace(0, deltat, N), rho_reducer_array)
# plt.xlabel('Time [s]')
# plt.ylabel('Density after reducer [kg/m^3]')
# plt.grid()

plt.subplot(133)
# plt.plot(np.linspace(0, deltat, N), p_press_array_ig/1e5, label='Ideal Gas')
plt.plot(np.linspace(0, deltat, N), p_press_array_td/1e5, label='Sys Ana Design')
plt.plot(np.linspace(0, deltat, N), p_press_array_td_w_ae/1e5, label='Sys Ana Design & Adiabatic Expansion')
plt.plot(np.linspace(0, deltat, N), p_press_array_jt/1e5, label='Real Gas / Joule Thompson')
plt.xlabel('Time [s]')
plt.ylabel('Pressurizer pressure [bar]')
plt.legend()
plt.grid()







#%% # 2 numerische integration von pressurizer tank entleerung und befüllung des propellant tanks
# mit der wesentlichen annahme, dass ideales gas, also T nach reducer = T vor reducer

# initialization

calc = False

if calc:
    ullage = 0.1
    vol_press = 0.041
    vol_press = 0.02868
    vol_proptank = 0.08*(1+ullage)
    vol_ullage = vol_proptank*ullage
    vol_prop = vol_proptank*(1-ullage)
    press = 'N2'
    rho_press_init = psi('D','P',p_init,'T',T_init,press)
    mass_press = rho_press_init*vol_press
    
    
    
    
    
    
    T_init = 293
    p_init = 300e5
    p_end = 80e5
    
    plot = False
    
    molar_mass = 28
    R = 8314
    
    p_reducer = 40e5
    
    m_dot_prop = 1
    T_prop = 293
    prop = 'C2H6O'
    rho_prop = psi('D','P',p_reducer,'T',T_prop,prop)
    
    # calculation
    
    vol_dot_prop = m_dot_prop/rho_prop
    
    N = 1000
    deltat = vol_prop / vol_dot_prop
    deltat = 120
    dt = deltat / N
    
    p_current = 300e5
    
    mass_press_array_ig = np.zeros(N)
    T_reducer_array_ig = np.zeros(N)
    rho_reducer_array_ig = np.zeros(N)
    p_press_array_ig = np.zeros(N)
    
    deltam = 0
    threshold = 6.04
    
    for ii in range(N):
        vol_dot_press = vol_dot_prop    
        # temperature pressurant after reducer equals temperature before reducer
        T_reducer = T_init*(p_current/p_init)**((kappa-1)/kappa)
        
        T_reducer = 201
        
        # T_reducer = T_end
        rho_reducer = p_reducer/T_reducer/(R/molar_mass)
        # rho_reducer = psi('D','P',p_reducer,'T',T_reducer,press) --> next step
        m_dot_reducer = vol_dot_press*rho_reducer
        mass_press -= m_dot_reducer*dt
        
        deltam += m_dot_reducer*dt
        # if deltam > threshold:
        #     break
        # if p_current < p_end:
        #     break
        
        rho_press = mass_press/vol_press
        p_current = p_init*(rho_press/rho_init)**kappa
        
        mass_press_array_ig[ii] = mass_press
        T_reducer_array_ig[ii] = T_reducer
        rho_reducer_array_ig[ii] = rho_reducer
        p_press_array_ig[ii] = p_current
        
    if plot:
    
        plt.figure(figsize=(12,8))
        
        plt.subplot(221)
        plt.plot(np.linspace(0, deltat, N), mass_press_array_ig)
        plt.xlabel('Time [s]')
        plt.ylabel('Mass Pressurizer [kg]')
        plt.grid()
        
        plt.subplot(222)
        plt.plot(np.linspace(0, deltat, N), T_reducer_array_ig)
        plt.xlabel('Time [s]')
        plt.ylabel('Temperature after reducer [K]')
        plt.grid()
        
        plt.subplot(223)
        plt.plot(np.linspace(0, deltat, N), rho_reducer_array_ig)
        plt.xlabel('Time [s]')
        plt.ylabel('Density after reducer [kg/m^3]')
        plt.grid()
        
        plt.subplot(224)
        plt.plot(np.linspace(0, deltat, N), p_press_array_ig/1e5)
        plt.xlabel('Time [s]')
        plt.ylabel('Pressurizer pressure [bar]')
        plt.grid()


# deutlich weniger pressurizer masse notwendig, als ich berechne






###################################################################
#%% altes zeug
###################################################################



#%% 3 ergäzung um joule thompson effekt durch coolprop
# plot = False

# mass_press = rho_press_init*vol_press

# mass_press_array_jt = np.zeros(N)
# T_reducer_array_jt = np.zeros(N)
# T_presstank_array_jt = np.zeros(N)
# p_press_array_jt = np.zeros(N)
# rho_presstank_array_jt = np.zeros(N)

# rho_press = rho_press_init
# p_current = p_init

# for ii in range(N):
#     vol_dot_press = vol_dot_prop    
#     # temperature pressurant after reducer equals temperature before reducer
#     T_presstank = psi('T','P',p_current,'D',rho_press,press)
#     h_presstank = psi('H','P',p_current,'D',rho_press,press)
#     T_reducer = psi('T','P',p_reducer,'H',h_presstank,press)
#     rho_reducer = psi('D','P',p_reducer,'T',T_reducer,press)
#     # rho_reducer = psi('D','P',p_reducer,'T',T_reducer,press) --> next step
#     m_dot_reducer = vol_dot_press*rho_reducer
#     mass_press -= m_dot_reducer*dt
    
#     rho_press = mass_press/vol_press
#     p_current = psi('P','T',T_presstank,'D',rho_press,press)
    
#     mass_press_array_jt[ii] = mass_press
#     T_reducer_array_jt[ii] = T_reducer
#     T_presstank_array_jt[ii] = T_presstank
#     p_press_array_jt[ii] = p_current
#     rho_presstank_array_jt[ii] = rho_press

# if plot:
    
#     plt.figure(figsize=(12,8))
    
#     plt.subplot(231)
#     plt.plot(np.linspace(0, deltat, N), mass_press_array_jt)
#     plt.xlabel('Time [s]')
#     plt.ylabel('Mass Pressurizer [kg]')
#     plt.grid()
    
#     plt.subplot(232)
#     plt.plot(np.linspace(0, deltat, N), T_reducer_array_jt)
#     plt.xlabel('Time [s]')
#     plt.ylabel('Temperature after reducer [K]')
#     plt.grid()
    
#     plt.subplot(233)
#     plt.plot(np.linspace(0, deltat, N), T_presstank_array_jt)
#     plt.xlabel('Time [s]')
#     plt.ylabel('Temperature pressurizer [K]')
#     plt.grid()
    
#     plt.subplot(234)
#     plt.plot(np.linspace(0, deltat, N), p_press_array_jt/1e5)
#     plt.xlabel('Time [s]')
#     plt.ylabel('Pressurizer pressure [bar]')
#     plt.grid()
    
#     plt.subplot(235)
#     plt.plot(np.linspace(0, deltat, N), rho_presstank_array_jt)
#     plt.xlabel('Time [s]')
#     plt.ylabel('Density Pressurizer [kg/m^3]')
#     plt.grid()
    
#     plt.subplot(236)
#     plt.plot(np.linspace(0, deltat, N), rho_presstank_array_jt/p_press_array_jt)
#     plt.xlabel('Time [s]')
#     plt.ylabel('Density Pressurizer [kg/m^3]')
#     plt.grid()
# #!!! offene fragen
# # ist es korrekt, dass die temperatur im pressurizer tank konstant bleibt???
# # oder kommt das durch eine falsche formulierung von mir?


# #%% 4 vergleich der ergebnisse

# plt.figure(figsize=(12,8))

# plt.subplot(221)
# plt.plot(np.linspace(0, deltat, N), mass_press_array, label='Ideal Gas')
# plt.plot(np.linspace(0, deltat, N), mass_press_array_jt, label='Real Gas CoolProp')
# plt.xlabel('Time [s]')
# plt.ylabel('Mass Pressurizer [kg]')
# plt.grid()
# plt.legend()

# plt.subplot(222)
# plt.plot(np.linspace(0, deltat, N), T_reducer_array, label='Ideal Gas')
# plt.plot(np.linspace(0, deltat, N), T_reducer_array_jt, label='Real Gas CoolProp')
# plt.xlabel('Time [s]')
# plt.ylabel('Temperature after reducer [K]')
# plt.grid()
# plt.legend()

# # plt.subplot(223)
# # plt.plot(np.linspace(0, deltat, N), T_press_array)
# # plt.xlabel('Time [s]')
# # plt.ylabel('Temperature pressurizer [K]')
# # plt.grid()

# plt.subplot(224)
# plt.plot(np.linspace(0, deltat, N), p_press_array/1e5, label='Ideal Gas')
# plt.plot(np.linspace(0, deltat, N), p_press_array_jt/1e5, label='Real Gas CoolProp')
# plt.xlabel('Time [s]')
# plt.ylabel('Pressurizer pressure [bar]')
# plt.grid()
# plt.legend()

# # 5 annahme überarbieten, dass T_2e = T_Ne --> s. power point progress update 05
# # stattdessen die temperature nach joule thompson berechnen und schaun was raus kommt

# # was ich stets vernachlässige bzw. annehme: adiabate tanks --> keine wärme fließt in die tanks