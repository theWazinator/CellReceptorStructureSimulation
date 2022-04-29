import math
import random
import numpy as np

def create_ligand_list(ligand_count, grid_point_length, grid_point_height):

    ligand_location_list = []

    for index in range(0, ligand_count):

        x_pos = random.randint(0, grid_point_length-1)
        y_pos = random.randint(0, grid_point_length-1)
        z_pos = random.randint(0, grid_point_height-1) # Location -1 is cell surface, location grid_point_length-1 is the layer before the ceiling

        location = [x_pos, y_pos, z_pos]

        ligand_location_list.append(location)

    return ligand_location_list

def ligand_in_receptor (ligand, receptor_list):

    for receptor in receptor_list:

        if ligand[0] == receptor[0] and ligand[1] == receptor[1]:

            return True

    return False

def get_equilibrium_current(ligands_received_time_list, time_steps):

    start_index = int(time_steps*0.5) # Take average of last 50% of current samples
    final_index = time_steps-1

    equilibrium_current = np.average(ligands_received_time_list[start_index:final_index])

    return equilibrium_current

def get_equilibrium_variance(ligands_received_time_list, time_steps):

    start_index = int(time_steps*0.5) # Take average of last 50% of current samples
    final_index = time_steps-1

    equilibrium_variance = np.var(ligands_received_time_list[start_index:final_index])

    return equilibrium_variance

# randomly-distributed receptors
def create_receptor_list_RANDOMDISTRIBUTION(grid_point_length, percent_coverage):

    # number_receptors = int(grid_point_length*grid_point_length*percent_coverage)

    receptor_location_list = []

    for x_index in range(0, grid_point_length):

        for y_index in range(0, grid_point_length):

            if random.random() < percent_coverage:

                location = [x_index, y_index]
                receptor_location_list.append(location)

    actual_percent_coverage = len(receptor_location_list) / (grid_point_length ** 2)

    print("Actual percent coverage: " +str(actual_percent_coverage))

    return receptor_location_list, actual_percent_coverage

def create_receptor_list_DISKDISTRIBUTION(grid_point_length, percent_coverage):

    N_disks = 20

    disk_radius = math.sqrt(grid_point_length*grid_point_length*percent_coverage/(N_disks*math.pi))

    print("Disk radius: " +str(disk_radius))

    receptor_location_list = []

    for index in range(0, N_disks):

        x_center = random.randint(0, grid_point_length-1)
        y_center = random.randint(0, grid_point_length-1)

        receptor_location_list.append([x_center, y_center])

        for x_coord in range(x_center-50, x_center+50):

            for y_coord in range(y_center-50, y_center+50):

                if math.sqrt((x_coord-x_center)**2 + (y_coord-y_center)**2) <= disk_radius:

                    receptor_location_list.append([x_coord % grid_point_length, y_coord % grid_point_length]) # wraparound the boundaries

    actual_percent_coverage = len(receptor_location_list) / (grid_point_length ** 2)

    print("Actual percent coverage: " + str(actual_percent_coverage))

    return receptor_location_list, actual_percent_coverage

def create_receptor_list_GRID(grid_point_length, percent_coverage):

    number_receptors = int(grid_point_length*grid_point_length*percent_coverage)

    receptors_per_line = math.sqrt(number_receptors)
    number_points_x = math.floor(receptors_per_line)
    number_points_y = math.ceil(receptors_per_line)
    spacing_x = int(grid_point_length/number_points_x)
    spacing_y = int(grid_point_length/number_points_y)

    spacing_x_list = []
    for index in range(0, grid_point_length):
        if index % spacing_x == 0 and len(spacing_x_list) < number_points_x:
            spacing_x_list.append(index)

    spacing_y_list = []
    for index in range(0, grid_point_length):
        if index % spacing_y == 0 and len(spacing_y_list) < number_points_y:
            spacing_y_list.append(index)

    print("Number Receptors: " +str(number_receptors))
    print("x-coordinates: " +str(spacing_x_list))
    print("y-coordinates: " +str(spacing_y_list))

    receptor_location_list = []

    for x_index in range(0, grid_point_length):

        if x_index in spacing_x_list:

            for y_index in range(0, grid_point_length):

                if y_index in spacing_y_list:

                    location = [x_index, y_index]
                    receptor_location_list.append(location)

    actual_percent_coverage = len(receptor_location_list) / (grid_point_length ** 2)

    print("Actual percent coverage: " + str(actual_percent_coverage))

    return receptor_location_list, actual_percent_coverage

def create_receptor_list_HEXAGON(grid_point_length, ratio):

    x_length_1 = 7*ratio
    x_length_2 = 3*ratio
    y_length = 3*ratio

    x_line_1_list = []

    x_length_1_boolean = True
    x_index = 0
    while x_index < grid_point_length:

        x_line_1_list.append(x_index)

        if x_length_1_boolean == True:
            x_index = x_index + x_length_1
            x_length_1_boolean = False

        else:
            x_index = x_index + x_length_2
            x_length_1_boolean = True

    x_line_2_list = []

    x_length_2_boolean = False
    x_index = 2
    while x_index < grid_point_length:

        x_line_2_list.append(x_index)

        if x_length_2_boolean == True:
            x_index = x_index + x_length_1
            x_length_2_boolean = False

        else:
            x_index = x_index + x_length_2
            x_length_2_boolean = True

    y_line_list_1 = []
    y_index = 0
    while y_index < grid_point_length:

        y_line_list_1.append(y_index)
        y_index = y_index+(2*y_length)

    y_line_list_2 = []
    y_index = 3
    while y_index < grid_point_length:
        y_line_list_2.append(y_index)
        y_index = y_index + (2*y_length)

    receptor_location_list = []

    for y_index in range(0, grid_point_length):

        if y_index in y_line_list_1:

            for x_index in x_line_1_list:

                location = [x_index, y_index]
                receptor_location_list.append(location)
                location = [x_index+1, y_index]
                receptor_location_list.append(location)
                location = [x_index, y_index+1]
                receptor_location_list.append(location)

        if y_index in y_line_list_2:

            for x_index in x_line_2_list:
                location = [x_index, y_index]
                receptor_location_list.append(location)
                location = [x_index + 1, y_index]
                receptor_location_list.append(location)
                location = [x_index, y_index + 1]
                receptor_location_list.append(location)

    actual_percent_coverage = len(receptor_location_list) / (grid_point_length ** 2)

    print("x_line_1_list: " +str(x_line_1_list))
    print("x_line_2_list: " + str(x_line_2_list))
    print("y_line_1_list: " + str(y_line_list_1))
    print("y_line_2_list: " + str(y_line_list_2))

    print("Actual percent coverage: " + str(actual_percent_coverage))

    return receptor_location_list, actual_percent_coverage










