"""
# Title: Self-assembly of polycubes (3D and 2D)

# Description: 
This script analyzes the self-assembly of 2D and 3D polyominoes. 
It processes input genotypes to return data on "valid" outputs, their frequencies, and complexities.
This code is broken down into modules and packages to support reusability.

[Author: Prarthana Agrawal]
[Date (last modified): 19 May 2025]
[Version: 1.1] Major fix over V0.1+

Older versions:
--------------------------------------------------------------------------------------------------------------------
    # 1.0 [created on 03 Mar 2025]
        - modification: when valid_shapes contained only 1 shape, the coordinated of the shape were being saved in different lines. To be consistent with multi-shape cases, I have changed the savetxt function to save all coordinates of one shape in one line.
        (minor): calculate complexity(species) and lz complexity only for valid shapes.
    
    # 0.4 [created on 04 Feb 2025]
        - errors: it used the wrong version of compare_polycubes (bbox was not a cube) leading to miscalculated valid solutions and their frequencies. I have fixed this by changing bbox to a cube. The error was spotted while plotting oeis upper bound on polycube number vs size.
        - modification:
            1. updated get_bbox function so compare_polycubes and valid_sol_checker are affected.
            2. simplify genotype rules
            3. calculate lz complexity and complexity (colors) within the main function.

    # 0.3 [created on 03 Feb 2025]
        - modification: 
            1. added unseeded assembly option. No direct change in main.py. 
            Changes reflect in assembly_func and valid_sol_checker.
            2. stopped storing all_orientations for 2D assembly.
        - any errors: 
            1. unseeded assembly does not have a path finder option
            2. if valid_shapes had only one entry [[(0,0,0)]] the savetxt option failed.

    # 0.2 [created on 28 Jan 2025]
        - modification: added exhaustive search for any number of tiles (even >2)
        - any errors: nsplit=0 was not defined for cased tot_splits < 1.

    # 0.1 [created on 13 Jan 2025] Built using 3D_V10+V22
        - modification: added another complexity measure (number of species)
        - any errors: None

--------------------------------------------------------------------------------------------------------------------
"""

#-------------------------------------- Import packages ----------------------------------------------------------#
# Inbuilt packages
import numpy as np
import random
import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.mplot3d import Axes3D
import itertools
import sys # required for command-line arguments

# User-defined packages
#from input import * # input values # use this when input.py is located in the same directory
from core import import_input, get_params, tile_identity_func # tile identities and system parameters
from core.validity_checks import valid_sol_checker # validity of solutions obtained in self-assembly
from utils import get_batch_range, save_all_genotypes, save_assembly_description, get_exh_genotypes
from symmetry import compare_polycubes
from plots import plot_all_cubes
from utils.genotype_utils import simplify_tiledict
from utils.binary_utils import tiledict_orientdict_to_binary, lempel_ziv_complexity

#-------------------------------------- Load input.py ----------------------------------------------------------#
# Get the input file path from the command-line argument
if len(sys.argv) < 2:
    print("Error: Please provide the path to input.py.")
    sys.exit(1)

input_file_path = sys.argv[1]
input_module = import_input(input_file_path)

# Extract the input values from the input.py file
n_tiles = input_module.n_tiles #int number of tile types / species
n_sides = input_module.n_sides #int number of side types / colors
neutral_sides = input_module.neutral_sides #list neutral sides
self_int_sides = input_module.self_int_sides #list self-interacting sides
n_rules = input_module.n_rules #int for random sampling or str 'all' for exhaustive search
dim = input_module.dim #int 2 or 3 for 2D or 3D assembly
assembly_type = input_module.assembly_type #str 'seeded' or 'unseeded'
Dmax = input_module.Dmax #int size cutoff for assembly
max_tiles = input_module.max_tiles #int size cutoff for assembly
kmax = input_module.kmax #int non-determinism check
tot_splits = input_module.tot_splits #int for 10 parallel splits , tot_splits = 10

# Check if the number of tiles is less than 26 (named A-Z)
if n_tiles > 26: raise NotImplementedError("n_tiles should be lesser than 26.")
if n_sides > 99: raise NotImplementedError("n_sides should be lesser than 99.")

# if the assembly is split into parallel batches, then get the split number
if tot_splits > 1: nsplit = int(sys.argv[2])
else: nsplit = 0 #! what does this mean? 
    
print(f"nsplit = {nsplit}, tot_splits = {tot_splits}")

# Random sampling or exhaustive search
if n_rules == 'all': # exhaustive search
    if dim != 2: raise NotImplementedError("Exhaustive search only possible for 2D assembly.")
    exhaustive_search = 'on'
    n_rules = (n_sides)**(2*dim*n_tiles) # total number of genotypes to assemble 
else: # random sampling
    exhaustive_search = 'off'
    
#------------------------------- Tile identities and system parameters ------------------------------------------#
params = get_params(n_tiles, n_sides, neutral_sides, self_int_sides, n_rules, dim)
assembly_settings = {
    'assembly_type': assembly_type,
    'Dmax': Dmax,
    'max_tiles': max_tiles,
    'kmax': kmax
}
parallel_run_config = {
    'tot_splits': tot_splits,
    'nsplit': nsplit #! is this nsplit even important?
}

#====================================================================================================================#
#----------------------------------------- Initialisation -----------------------------------------------------------#
#====================================================================================================================#

global UBD, ND, valid
UBD = 0; ND = 0; valid = 0 # number of unbounded (UBD), non deterministic (ND) and valid solutions
sol_stats = dict({'UBD': UBD, 'ND': ND, 'valid': valid}) # solution statistics

# -------------------------------------- Lists to store data -------------------------------------------------------#
frequency = list() # list of frequency of each valid polycube
complexity_list = list() # list of complexity of each valid polycube == min number of colors used
complexity_species_list = list() # alternate complexity measure == min number of species used
lz_complexity_list = list() # list of Lempel-Ziv complexity of each valid polycube

#all_shapes = list() # list of all polycubes generated
valid_shapes = list() # list of all "valid" polycubes
shape_type = list() # list of shape "types" (numbered in order of appearance) generated by each genotype
all_genotypes = list() # list of all genotypes
all_orientations = list() # list of all orientations of tiles
shape_type_dict = dict() # identify a shape to assign a shape "type" to it

# randomly assign a master seed
master_seed = random.randint(0, n_rules) # take a random seed value
random.seed(master_seed) #! To regenerate results -- doesn't work for unseeded assembly

#====================================================================================================================#
#---------------------------------- Get range of nrules to assemble in this run -------------------------------------#
#====================================================================================================================#
# if parallel batches are created, find the indices for this particular nsplit
indices, file_name = get_batch_range(params, parallel_run_config)
start_index, end_index = indices
batch_size = round((end_index-start_index)/ 10) # round off to nearest-10

# if exhaustive search is on, get list of all genotypes in order
if exhaustive_search == 'on':
    if dim != 2: raise NotImplementedError("Exhaustive search only possible for 2D assembly.")

    # Generate valid genotype combinations per tile
    side_types = params['side_types']
    genotype_combinations_per_tile = [list(side) + ['00', '00'] for side in itertools.product(side_types, repeat=4)]
    num_geno_per_tile = len(genotype_combinations_per_tile)

    # Generate only the required slice of genotype combinations
    indices_range_for_this_nsplit = list(itertools.islice(
        itertools.product(range(num_geno_per_tile), repeat=n_tiles),
        start_index, end_index))
    
    kr = 0  # Counter for pair_split
    
#====================================================================================================================#
#------------------------- Assemble polycubes for each genotype in this selected range ------------------------------#
#====================================================================================================================#

for rule_num in range(start_index, end_index):
    if rule_num % batch_size == 0: print("#", rule_num)

    #==============================================================================================================#
    #---------------------------------------- Get tile_dict, orient_dict ------------------------------------------#
    #==============================================================================================================#
    
    # --------- Random sampling to generate a tile and orientation dictionary -------#
    if exhaustive_search == 'off':
        tile_dict, orient_dict = tile_identity_func(params)

    # --------- Exhaustive search to generate a tile and orientation dictionary -------#
    elif exhaustive_search == 'on':

        tile_dict = {}
        orient_dict = {}

        tile_labels = (chr(i) for i in range(65, 65 + n_tiles))  # Assign labels A, B, C, ...
        index_tuple = indices_range_for_this_nsplit[kr]

        for i, tile in enumerate(tile_labels):
            tile_dict[tile] = genotype_combinations_per_tile[index_tuple[i]]
            orient_dict[tile] = "OOOOUU"  # For 2D, orientations don't matter
            if dim == 3: raise NotImplementedError("orientation not implemented.")

        kr += 1  # Move to the next genotype set

    else: raise ValueError("Indicate whether exhaustive search or random sampling.")
    #print(tile_dict, "\n", orient_dict)

    #---------- Manually check output for a tile_dict and orient_dict ------------#
    #tile_dict = dict({'A': ['00', '03', '02', '01', '03', '02'], 'B': ['00', '01', '00', '03', '00', '02']} )
    #orient_dict = dict({'A': 'ODLDUD', 'B': 'RURURL'})
    
    # store the genotype and orientations (store original, non-simplified versions)
    all_genotypes.append(list(tile_dict.values()))
    if dim == 3: all_orientations.append(list(orient_dict.values()))

    # Simplify tile_dict before assembling
    # 1. '00' all neutral sides 2. '00' sides with no partners 3. renumber sides
    tile_dict = simplify_tiledict(params, tile_dict)

    #==============================================================================================================#
    #------------------------------------- Main assembly for this genotype ----------------------------------------#
    #==============================================================================================================#
    
    # check if the resulting phenotype is valid or not
    output, tile_coord, picked_tiles, complexity, sol_stats = valid_sol_checker(params, assembly_settings, tile_dict, orient_dict, sol_stats)

    # comp(species) and lz-complexity calculation
    if output == 1: #? new addition --  no point calculating complexity for invalid shapes
        complexity_species = len(set(picked_tiles)) # complexity measure: min number of tile species required

        # Remove tiles with all sides as '00' before calculating lz complexity
        nonzero_tile_dict = {k: v for k, v in tile_dict.items() if not all(side == '00' for side in v)}
        nonzero_orient_dict = {k: orient_dict[k] for k in nonzero_tile_dict}
        binary_genotype = tiledict_orientdict_to_binary(n_sides, nonzero_tile_dict, nonzero_orient_dict)
        lz_complexity = lempel_ziv_complexity(binary_genotype)
            
    #==============================================================================================================#
    #----------------------- Find out if this is a new shape, store freq and complexities--------------------------#
    #==============================================================================================================#
    
    # default value: the new phenotype is not same as any old encountered phenotype
    match = 'no'

    # plot all shapes produced valid or not
    #%matplotlib widget
    #fig = plt.figure()
    #ax1 = fig.add_subplot(111, projection='3d')
    #ax1.view_init(elev=-149, azim=138)
    #plot_all_cubes(params, ax1, tile_coord, picked_tiles, cube_outline='False', axes_lines='False')
    
    if output == 1: # shape is valid
        # plot all valid shapes
        #fig = plt.figure()
        #ax1 = fig.add_subplot(111, projection='3d')
        #ax1.view_init(elev=-149, azim=138)
        #plot_all_cubes(params, ax1, tile_coord, picked_tiles, cube_outline='False', axes_lines='False')

        if len(valid_shapes) > 0: # second entry onwards
            for i in range(len(valid_shapes)):
                output1 = compare_polycubes(tile_coord, valid_shapes[i]) 
                if output1 == 1: 
                    frequency[i] = frequency[i] + 1
                    complexity_list[i] = min(complexity_list[i], complexity)
                    complexity_species_list[i] = min(complexity_species_list[i], complexity_species)
                    lz_complexity_list[i] = min(lz_complexity_list[i], lz_complexity)
                    shape_type.append(int(shape_type_dict[np.array(valid_shapes[i]).tobytes()]))
                    match = 'yes'
                    break

            if match == 'no':
                valid_shapes.append(tile_coord)
                frequency.append(1)
                complexity_list.append(complexity)
                complexity_species_list.append(complexity_species)
                lz_complexity_list.append(lz_complexity)
                shape_type.append(max(shape_type)+1)
                shape_type_dict[np.array(tile_coord).tobytes()] = max(list(shape_type_dict.values()))+1

        else: # for the first entry
            valid_shapes.append(tile_coord)
            frequency.append(1)
            complexity_list.append(complexity)
            complexity_species_list.append(complexity_species)
            lz_complexity_list.append(lz_complexity)
            shape_type_dict[np.array(tile_coord).tobytes()] = 1
            shape_type.append(int(1))
    else:
        #all_shapes.append(np.array([0])) # 0 implies that the shape is invalid
        shape_type.append(int(0))

plt.show()

# ------------------------------------- Save data files ----------------------------------------------------------#

# Modify file names to include the directory path

#debian
#path_name = "/media/agrawalp/221ceb7e-aa6b-4034-b2ee-bb33de6397d3/polyominoes_new/polycube/runs/{}D/".format(dim)

#if exhaustive_search == 'on':
#    folder_name = f"{n_tiles}s{n_sides}c_exh/data_files/"
#else:
#    folder_name = f"{n_tiles}s{n_sides}c/nrules_10^{int(np.log10(n_rules))}/data_files/" 
    
#file_name = path_name + folder_name + file_name 

#cluster
file_name = "data_files/" + file_name #windows

if dim == 3: save_orientations = True 
else: save_orientations = False

save_all_genotypes(file_name, all_genotypes, all_orientations, save_orientations)
save_assembly_description(file_name, params, assembly_settings, parallel_run_config, master_seed, sol_stats)

np.savetxt(file_name + 'shape_type.txt', shape_type, fmt='%d')
np.savetxt(file_name + 'frequency.txt', frequency, fmt='%d')
np.savetxt(file_name + 'complexity.txt', complexity_list, fmt='%d')
np.savetxt(file_name + 'complexity_species.txt', complexity_species_list, fmt='%d')
np.savetxt(file_name + 'lz_complexity.txt', lz_complexity_list, fmt='%d')

if len(valid_shapes) == 1: # if valid_shapes = [[(0,0,0)]], then savetxt interprets it as a 3D array
    with open(file_name + 'valid_shapes.txt', 'w') as f:
        f.write(str(valid_shapes[0]) + '\n')
else:
    valid_shapes_to_save = np.array(valid_shapes, dtype=object)
    np.savetxt(file_name + 'valid_shapes.txt', valid_shapes_to_save, fmt='%s')