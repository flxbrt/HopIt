# -*- coding: utf-8 -*-
"""
Title: main.py
Project: HopIt - Hopper Iterative System Analysis Tool
Author: @flxbrt
Version: 2.1

Description:
    Entry point for iterative hopper system analysis. Loads configuration,
    initializes combustion parameters, and runs a convergence loop to estimate
    total system mass through subsystem evaluations.

Usage:
    python main.py

Dependencies:
    - core.functions
    - core.casadi_core
"""

import os
import time
import core.functions as fc
from core.casadi_core import solve_tanks

# Get the directory of the current script (for relative paths)
path = os.path.dirname(__file__)

# System configuration parameters (could also be loaded from a JSON or YAML file)
config = {
    'Flight Time [s]': 80,
    'Oxidizer Pump': False,
    'Fuel Pump': False,
    'Oxidizer': 'O2',
    'Fuel': 'C2H5OH',
    'Fuel Coolprop': 'C2H6O',
    'Fuel Self Pressurised': False,
    'Oxidizer Self Pressurised': False
}

# Control flags for printing and storing the final system
sys_print = False
store_system = False

# Load initial system state based on the config
system = fc.load_system(config, path)

# Set up combustion properties (e.g. equilibrium, fuel mixture)
comb = fc.set_combustion()

#%% Main iterative loop

# Compute ullage (gas volume margin) for pressurization
system = fc.calculate_ullage(system)

counter = 0
system_mass_prev = 1  # Initial value to trigger the while-loop
system_mass = fc.get_system_mass(system)

# Estimate exhaust velocity once before loop (isentropic, max conditions)
system = fc.calc_isentropic_exhaust_velocity(system, comb, 'max')

threshold = 0.1  # [kg] — convergence threshold

t_start = time.time()

# Iterative loop: repeat calculations until mass converges
while abs(system_mass_prev - system_mass) > threshold:
    system_mass_prev = system_mass

    # 1) Estimate propellant mass using Ziolkowsky equation
    system = fc.ziolkowsky(system)

    # 2) Estimate tank mass based on propellant volume
    system = solve_tanks(system)

    # 3) Evaluate thrust chamber (mass, thrust, specific impulse)
    system = fc.solve_thrust_chamber(system, comb)

    # 4) Calculate battery and pump sizing (if applicable)
    system = fc.compute_pumps_and_batteries(system)

    # 5) Calculate structural mass (requires max thrust info)
    system = fc.solve_structure_frame(system)

    # 6) Add a dry mass margin (for robustness)
    system = fc.add_mass_margin(system)

    # Update system mass and print iteration info
    system_mass = fc.get_system_mass(system)
    counter += 1
    print(f'Iteration {counter} done')
    print(f'Relative deviation: {(system_mass - system_mass_prev) / system_mass_prev * 100:.3f} [%]\n')

# Compute valve sizing once mass has converged
system = fc.compute_valves(system)

t_end = time.time()

print(f'Converged in {t_end - t_start:.2f} [s] and {counter} Iterations')

#%% Optional system output
if sys_print:
    fc.print_system(system)

#%% Optional system storage
if store_system:
    fc.write_system(system)