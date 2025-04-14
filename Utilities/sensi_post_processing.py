# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 07:31:34 2024

@author: felix
"""



import os
import json
import numpy as np
import matplotlib.pyplot as plt



def load_system(direc, file_name):
    # config_direc = 'C:/Users/felix/Desktop/System Analysis/Current/'
    with open(direc + file_name, 'r') as system_file:
        system = system_file.read()
    system = json.loads(system)
    return system

def get_system_mass(system, ignore=list()):
    mass = 0
    for subsystem in system.keys():
        if subsystem != 'Config':
            for component in system[subsystem].keys():
                if component in ignore:
                    continue
                mass += system[subsystem][component]['Mass']['Value']
    return mass

def get_sensitivity(x_values, y_values):
    
    index_max = y_values.index(max(y_values))
    index_min = y_values.index(min(y_values))
    
    return (y_values[index_max] - y_values[index_min])/(x_values[index_max] - x_values[index_min])



direc = 'C:/Users/felix/Desktop/PDR/Sensitivity_analysis/'
save_direc = direc + 'plots/'


# TODO
key = 'LOx'

lim = True
val_sub = 0

if key == 'Flighttime':
    index = 3
    x_axis_label = 'Flight Time [s]'
    unit = '[kg/s]'
    lim = False
elif key == 'ROF':
    index = 2
    x_axis_label = 'ROF [-]'
    unit = '[kg/-]'
elif key == 'dry_mass':
    index = 3
    x_axis_label = 'Dry Mass Deviation [kg]'
    unit = '[kg/kg]'
elif key == 'thrust_chamber_efficiency':
    index = 3
    x_axis_label = 'Thrust Chamber Efficiency [%]'
    unit = '[kg/%]'
elif key == 'pcc':
    index = 1
    x_axis_label = 'Chamber Pressure @ 100% Thrust [bar]'
    unit = '[kg/bar]'
elif key == 'fueltank':
    index = 2
    x_axis_label = 'Fuel Line Pressure Loss [bar]'
    unit = '[kg/bar]'    
    val_sub = 20
elif key == 'oxtank':
    index = 1
    x_axis_label = 'Ox Line Pressure Loss [bar]'
    unit = '[kg/bar]'    
    val_sub = 20
elif key == 'pressurant_end':
    index = 3
    x_axis_label = 'Pressurant End Pressure [bar]'
    unit = '[kg/bar]'    
elif key == 'pressurant_init':
    index = 3
    x_axis_label = 'Pressurant Init Pressure [bar]'
    unit = '[kg/bar]'    
elif key == 'LOx':
    index = 1
    x_axis_label = 'LOx temperature [K]'
    unit = '[kg/K]'    
    
    
    

    








all_files = os.listdir(direc)



files = list()

for file_name in all_files:
    if key in file_name:
        files.append(file_name)
    

results = list()
x_values = list()



for file_name in files:
    system = load_system(direc, file_name)
    results.append(get_system_mass(system))
    
    
    
    
# TODO
    # ROF
    # val = float(file_name.split('_')[-1].split('j')[0][:-1])
    # flight time
    # val = float(file_name.split('_')[-1].split('s')[0])
    # dry mass & thrust chamber efficiency
    # val = float(file_name.split('_')[-1].split('.')[0])
    # pcc / fueltank / oxtank / pressurant end & init
    # val = float(file_name.split('_')[-1].split('b')[0])
    # LOx temp
    val = float(file_name.split('_')[-1].split('K')[0])

    
    
    x_values.append(val)     #.append(get_system_mass(system, ['Pressurant', 'Fuel', 'Oxidizer', 'Pressurant Tank', 'Oxidizer Tank', 'Fuel Tank', 'Frame', 'Margin']))
 
# TODO
results = [r for _,r in sorted(zip(x_values, results))] 
x_values = sorted(x_values)
# x_values = [x for _,x in sorted(zip(results,x_values))]
# results = sorted(results)

#%%






x_values = [x - val_sub for x in x_values]

# font = {'size': 30}

# Create a figure and an axes
fig, ax = plt.subplots(figsize=(10, 6))  # Optional: Specify the size of the figure

# Plot the data
ax.plot(x_values, results, marker='o', linestyle='-')#, color='b')#, label='Population')

# index = 3
ax.plot(x_values[index], results[index], 'ro', marker='D', markersize=10)

# ax.font(font)

# Add labels and title

ax.set_xlabel(x_axis_label, fontsize = 16)
ax.set_ylabel('System Mass [kg]', fontsize = 16)
sensi = get_sensitivity(x_values, results)
ax.set_title(f'Sensitivity is {sensi:.2f} {unit}', fontsize = 16)

ax.tick_params(axis='both', which='major', labelsize=14)


if lim == True:
    ax.set_ylim(165, 290)

# Add a legend
# ax.legend(fontsize = 30)

# Add a grid
ax.grid(True)

# Optionally, add annotations
# for i, (year, pop) in enumerate(zip(years, population)):
    # ax.annotate(f'{pop}', (year, pop), textcoords="offset points", xytext=(0,10), ha='center')

# Show the plot
plt.show()
plt.savefig(save_direc + key + '.svg')