# assumption about how molecule is instantly transported into cell and how each transmembrane protein acts entirely independently (i.e. there is no signalling and coordinated opening/closing of channels or movement towards/away from each other.
# assume no collisions between molecules
# mention how each receptor takes up a 5x5 nm area, since receptors are fixed-size in real life, so we keep this constant
# Equilibrium concentration after burn-in time

import math
import random
import pandas as pd
import numpy as np

from helper_methods import *

# Model parameters
ligand_diameter = 1e-9 # m
k_B = 1.38e-23 # J/K
nu = 1e-3 # kg/(m*s)
temperature = 300 # K

A_n = 1e23*6.022

diffusion = (k_B*temperature)/(6*math.pi*nu*ligand_diameter) # m^2/s
delta_x = 5*ligand_diameter # the distance between each grid point
delta_t = (3*math.pi*nu*math.pow(delta_x, 3))/(k_B*temperature) # time to move one body radius
grid_point_length = 300
grid_point_height = 50
base_time_steps = 500
grid_point_volume = grid_point_length*grid_point_length*grid_point_height
volume = grid_point_volume*math.pow(delta_x, 3)

random_seed = 12345678

random.seed(random_seed)

print("Ligand Size: " +str(ligand_diameter))
print("Diffusion Constant: " +str(diffusion))
print("Base Time steps: " + str(base_time_steps))
print("Delta_t: " +str(delta_t))

print("Grid point length count: " +str(grid_point_length))
print("Grid point area count: " +str(grid_point_length*grid_point_length))
print("Grid point volume count: " +str(grid_point_volume))
print("Grid space volume: " +str(volume))
print("Maximum concentration limit: " +str(1/math.pow(delta_x, 3))+ " ligands/m^3")
print("Maximum concentration limit: " +str(1/(math.pow(delta_x, 3)*A_n))+ " M")

mole_list = [1, 0.8, 0.6, 0.4, 0.2, 0.01, 0.008, 0.006, 0.004, 0.002]
percent_coverage_list = [0.01, 0.005, 0.0025, 0.00125, 0.000625]

number_ligands_to_simulate_list = []

for concentration in mole_list:

    num_molecules_liter = concentration*A_n
    num_molecules_volume = num_molecules_liter*volume
    number_ligands_to_simulate_list.append(int(num_molecules_volume))

current_dict = {}
variance_dict = {}

for percent_coverage in percent_coverage_list:

    current_list = []
    variance_list = []

    for molarity_index in range(0, len(mole_list)):

        molarity = mole_list[molarity_index]
        time_steps = int(base_time_steps/molarity)

        print("Molarity: "+str(molarity))
        print("Area coverage: " +str(percent_coverage))
        print("Number of ligands to simulate: " +str(number_ligands_to_simulate_list[molarity_index]))
        print("Number of time steps: " +str(time_steps))
        print("Total simulation time: " + str(time_steps * delta_t))

        # Create ligand list
        ligand_list = create_ligand_list(number_ligands_to_simulate_list[molarity_index], grid_point_length, grid_point_height)

        # Create receptor location list
        receptor_list, actual_percent_coverage = create_receptor_list_RANDOMDISTRIBUTION(grid_point_length, percent_coverage)

        # Only need to check if in receptor if z value is 0

        print("Begin Simulation")

        ligands_received_time_list = []

        for t in range(0, time_steps):

            if t%(int(time_steps/10)) == 0:

                print("Current timestep: " +str(t))

                if t > 0:
                    print("Average of last "+ str(int(time_steps/10)) +": " +str(np.average(ligands_received_time_list[-int(time_steps/10):-1])))

            ligands_received = 0

            for ligand in ligand_list:

                rand_dir = random.randint(0, 5)

                if rand_dir == 0:

                    ligand[2] = ligand[2] - 1 # go down

                elif rand_dir == 1:

                    ligand[2] = ligand[2] + 1 # go up

                elif rand_dir == 2:

                    ligand[1] = ligand[1] - 1 # go backward

                elif rand_dir == 3:

                    ligand[1] = ligand[1] + 1 # go forward

                elif rand_dir == 4:

                    ligand[0] = ligand[0] - 1 # go left

                elif rand_dir == 5:

                    ligand[0] = ligand[0] + 1 # go right

                if ligand[0] < 0:

                    ligand[0] = grid_point_length - 1 # periodic boundary

                if ligand[0] >= grid_point_length:

                    ligand[0] = 0 # periodic boundary

                if ligand[1] < 0:

                    ligand[1] = grid_point_length - 1 # periodic boundary

                if ligand[1] >= grid_point_length:

                    ligand[1] = 0 # periodic boundary

                if ligand[2] >= grid_point_height:

                    ligand[0] = random.randint(0, grid_point_length-1) # Place the ligand at a random location at the ceiling
                    ligand[1] = random.randint(0, grid_point_length-1)
                    ligand[2] = grid_point_height - 1

                if ligand[2] < 0:

                    if ligand_in_receptor(ligand, receptor_list) == True:

                        ligands_received = ligands_received + 1

                        ligand[0] = random.randint(0, grid_point_length - 1)  # Place the ligand at a random location at the ceiling
                        ligand[1] = random.randint(0, grid_point_length - 1)  # The total number of ligands (and concentration) is maintained
                        ligand[2] = grid_point_height - 1

                    else:

                        ligand[2] = 0  # ligand "bounces off" cell surface

            ligands_received_time_list.append(ligands_received)

        equilibirum_current = get_equilibrium_current(ligands_received_time_list, time_steps)

        equilibrium_variance = get_equilibrium_variance(ligands_received_time_list, time_steps)

        print("For molarity " +str(molarity)+ " M and area coverage" +str(actual_percent_coverage)+ "%: " +str(equilibirum_current) +"ligands absorbed per " +str(delta_t)+ "s and " +str(equilibrium_variance)+ " variance.", flush=True)

        current_list.append(equilibirum_current)
        variance_list.append(equilibrium_variance)

    current_dict[str(actual_percent_coverage)] = current_list
    variance_dict[str(actual_percent_coverage)] = variance_list

print("End Simulation")

current_df = pd.DataFrame.from_dict(current_dict)
current_df.to_csv(path_or_buf="RANDOMDISTRIBUTION_Current.csv")

variance_df = pd.DataFrame.from_dict(variance_dict)
variance_df.to_csv(path_or_buf="RANDOMDISTRIBUTION_Variance.csv")