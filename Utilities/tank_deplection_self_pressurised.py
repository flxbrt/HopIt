# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 11:58:27 2024

@author: felix
"""



import time
import numpy as np
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI as psi


def rk4_e(f, y, h, t, *args):
    # runge kutte 4th order explicit
    tk_05 = t + 0.5*h
    yk_025 = y + 0.5 * h * f(t, y, *args)
    yk_05 = y + 0.5 * h * f(tk_05, yk_025, *args)
    yk_075 = y + h * f(tk_05, yk_05, *args)
    
    return y + h/6 * (f(t, y, *args) + 2 * f(tk_05, yk_025, *args) + 2 * f(tk_05, yk_05, *args) + f(t+h, yk_075, *args))


def get_rho(T, Q):
    return psi('D', 'T', T, 'Q', Q, 'N2O')    

def get_p(T):
    # no difference, whether Q = 0 or Q = 1
    return psi('P', 'T', T, 'Q', 0, 'N2O')

def get_u(T, Q):
    return psi('U', 'T', T, 'Q', Q, 'N2O')
    
def get_h(T, Q):
    return psi('H', 'T', T, 'Q', Q, 'N2O')
    
def ode(t, y, *args):
    
    m_out = args[0]
    h_out = args[1] # saturated liquid according to paper p.8 
    Q_in_l = 0 #args[2] # from air into the wall
    Q_out_l = 0 #args[3] # from wall / out of wall into the fluid inside the tank
    Q_in_g = 0 #args[4]
    Q_out_g = 0 #args[5]
    m_w = 0 #args[6]
    c_w = 0 #args[7]
    
    # p1_f = y[0]             # Pressure in pipe 1
    
    
    dm_tot = -m_out
    dU_tot = -m_out*h_out + Q_out_l + Q_out_g
    '''
    ist m_w die gesamtmasse der wand? oder nur der jeweilig benetzte teilabschnitt?
    '''
    dTw_l = 0 #(Q_in_l - Q_out_l) / m_w/c_w
    dTw_g = 0 #(Q_in_g - Q_out_g) / m_w/c_w
    
    return np.array([dm_tot, dU_tot, dTw_l, dTw_g])

'''
def convective_heat():
    Ra = 0
    Nu = 0
    h = 0
    Q_dot = h*A*deltaT
'''
    
def find_temperature(V_tank, T_prev, U_tot, m_tot):
    # x = 0.5 # kenn ich x oder ergibt sich das aus dem gleichungssystem? --> 2 gleichungen und zwei unbekannte?!
    
    err = 1
    T = T_prev
    
    while err > 1e-4: # 1e-3 corresponds to 1l inaccuracy
        u_l = get_u(T, 0)
        u_g = get_u(T, 1)
        rho_l = get_rho(T, 0)
        rho_g = get_rho(T, 1)
        x = (U_tot/m_tot - u_l) / (u_g - u_l)
        V_tank_calc = m_tot*((1-x)/rho_l + x/rho_g)
        if V_tank_calc - V_tank > 0:
            sign = 1
        elif V_tank_calc - V_tank < 0:
            sign = -1
        else:
            sign = 0
        err = abs(V_tank_calc - V_tank) # error computation
        T = T + sign*0.01 # correction
    
    return T, x
       
def find_volume_fraction(x, T, m_tot, V_tot):
    rho_g = get_rho(T, 1)
    w = x*m_tot/rho_g/V_tot
    return w


def main_loop(y, T, timesteps, V_tank):
    
    sol = np.zeros((6, len(timesteps)+1))
    sol[:, 0] = np.array([y[0], y[1], y[2], y[3], T, get_p(T)])
    
    fractions = np.zeros((2, len(timesteps)))
    
    t_start = time.time()
    for t_index, t in enumerate(timesteps):
        # print(t_index)
        
        T, x = find_temperature(V_tank, T, y[1], y[0])
        
        if x>=1:
            break
        
        w = find_volume_fraction(x, T, y[0], V_tank)
        
        fractions[:, t_index] = np.array([x, w])
        
        # Q = 0 gilt nur solange Flüssigkeit aus dem Tank fließt
        h_out = get_h(T, 0)
        
        # Q_in_l = 0
        # Q_out_l = 0
        # Q_in_g = 0
        # Q_out_g = 0
        # if heat_ext:
        #     pass

        args = [m_out, h_out]#, Q_in_l, Q_out_l, Q_in_g, Q_out_g, m_w, c_w]
        
        y = rk4_e(ode, y, h, t, *args)
        
        sol[:4, t_index+1] = y
        sol[4:, t_index+1] = np.array([T, get_p(T)])
        
    t_end = time.time()
    print(f'Numerical Integration took {(t_end-t_start):.2f} [s]')
    
    return sol, fractions, t_index


def plot_transient(timesteps, t_index, sol, fractions):
    
    plt.figure(figsize=(16,8))
    plt.subplot(231)
    plt.plot(timesteps[:t_index], sol[0,:t_index], label='Total Mass')
    plt.plot(timesteps[:t_index], sol[0,:t_index]*(1-fractions[0,:t_index]), label='Liquid Mass')
    plt.plot(timesteps[:t_index], sol[0,:t_index]*fractions[0,:t_index], label='Gas Mass')
    plt.xlabel('Time [s]')
    plt.ylabel('Mass [kg]')
    plt.grid()
    plt.legend()

    plt.subplot(232)
    plt.plot(timesteps[:t_index], fractions[0,:t_index], label='Mass')
    plt.plot(timesteps[:t_index], fractions[1,:t_index], label='Volume')
    plt.xlabel('Time [s]')
    plt.ylabel('Gas Fraction [-]')
    plt.grid()
    plt.legend()

    plt.subplot(234)
    plt.plot(timesteps[:t_index], sol[4,:t_index])
    plt.xlabel('Time [s]')
    plt.ylabel('Temperature [K]')
    plt.grid()

    plt.subplot(235)
    plt.plot(timesteps[:t_index], sol[5,:t_index]/1e5)
    plt.xlabel('Time [s]')
    plt.ylabel('Pressure [bar]')
    plt.grid()

    plt.subplot(236)
    plt.plot(timesteps[:t_index], get_rho(sol[4,:t_index], 0))
    plt.xlabel('Time [s]')
    plt.ylabel('Density [kg/m^3]')
    plt.grid()


def plot_gas_fraction(sol, vol_frac, indices):
    
    pressure = np.zeros(len(vol_frac))
    for num_index, t_index in enumerate(indices):
        pressure[num_index] = sol[5, int(t_index), num_index]
    
    plt.plot(vol_frac, pressure/1e5)
    plt.xlabel('Initial Volumetric Gas Fraction')
    plt.ylabel('End Pressure [bar]')
    plt.grid()



exec_ = 3



#%% 1) transient tank depletion
if exec_ == 1:
    V_tank = 10e-3
    T0 = 293
    m_out = 1
    
    w0 = 0.1 # volumetric amount of vapor
    
    rho_l0 = get_rho(T0, 0)
    m_l0 = (1-w0)*V_tank*rho_l0
    rho_g0 = get_rho(T0, 1)
    m_g0 = w0*V_tank*rho_g0
    m_tot0 = m_l0 + m_g0
    
    u_l0 = get_u(T0, 0)
    u_g0 = get_u(T0, 1)
    U_tot0 = u_l0 * m_l0 + u_g0 * m_g0
    
    Tw_l0 = 0
    Tw_g0 = 0
    
    h = 1e-2
    duration = 100
    N = int(duration/h)
    timesteps = np.linspace(0, duration, N+1)
    
    y0 = np.array([m_tot0, U_tot0, Tw_l0, Tw_g0])
    
    y = y0
    T = T0
    
    # sol = np.zeros((6, len(timesteps)+1))
    # sol[:, 0] = np.array([y0[0], y0[1], y[2], y[3], T, get_p(T)])
    # fractions = np.zeros((2, len(timesteps)))
    
    sol, fractions, t_index = main_loop(y, T, timesteps, V_tank)
    
    plot_transient(timesteps, t_index, sol, fractions)

#%% 2) impact of initial volumetric gas fraction on end pressure
if exec_ == 2:
    V_tank = 10e-3
    T = 298
    m_out = 1
    
    Tw_l0 = 0
    Tw_g0 = 0
    
    rho_l0 = get_rho(T, 0)
    rho_g0 = get_rho(T, 1)
    u_l0 = get_u(T, 0)
    u_g0 = get_u(T, 1)
    
    h = 1e-2
    duration = 100
    N = int(duration/h)
    timesteps = np.linspace(0, duration, N+1)
    
    vol_frac = np.linspace(0.01, 0.8, 10)
    
    sol = np.zeros((6, len(timesteps)+1, len(vol_frac)))
    fractions = np.zeros((2, len(timesteps), len(vol_frac)))
    indices = np.zeros(len(vol_frac))
    
    for w_index, w in enumerate(vol_frac):
        
        m_l0 = (1-w)*V_tank*rho_l0
        m_g0 = w*V_tank*rho_g0
        m_tot0 = m_l0 + m_g0
        
        U_tot0 = u_l0 * m_l0 + u_g0 * m_g0
        
        y = np.array([m_tot0, U_tot0, Tw_l0, Tw_g0])
        
        sol_, fractions_, t_index = main_loop(y, T, timesteps, V_tank)
        
        sol[:, :, w_index] = sol_
        fractions[:, :, w_index] = fractions_
        indices[w_index] = t_index

    plot_gas_fraction(sol, vol_frac, indices)
    
#%% 3) identifying initial temperature - gas mas fraction pairs to end up with a certain end pressure
if exec_ == 3:
    
    V_tank = 10e-3
    m_out = 1
    
    Tw_l0 = 0
    Tw_g0 = 0
    
    h = 1e-1
    duration = 100
    N = int(duration/h)
    timesteps = np.linspace(0, duration, N+1)
    
    p_end_threshold = 35e5
    temperature = np.linspace(290, 300, 10)
    
    vol_frac = np.zeros(len(temperature))
    
    for index_T, T in enumerate(temperature):
        
        w = 0
        
        rho_l0 = get_rho(T, 0)
        rho_g0 = get_rho(T, 1)
        u_l0 = get_u(T, 0)
        u_g0 = get_u(T, 1)
        
        p_end = 0
        
        while p_end<p_end_threshold and w<=0.8:            
            w += 0.01
            print(f'Evalating p_end for outer loop {index_T+1} of {len(temperature)} at {w=}')
            m_l0 = (1-w)*V_tank*rho_l0
            m_g0 = w*V_tank*rho_g0
            m_tot0 = m_l0 + m_g0
            
            U_tot0 = u_l0 * m_l0 + u_g0 * m_g0
            
            y = np.array([m_tot0, U_tot0, Tw_l0, Tw_g0])
            
            sol_, fractions_, t_index = main_loop(y, T, timesteps, V_tank)
            
            p_end = sol_[5, t_index]
            
        vol_frac[index_T] = w

ullage = vol_frac/(1-vol_frac)

plt.plot(temperature, vol_frac, label='Initial Gas Volume')
plt.plot(temperature, ullage, label='Ullage for Sys Ana')
plt.xlabel('Initial Temperature')
plt.ylabel('Fraction')
plt.grid()
plt.legend()