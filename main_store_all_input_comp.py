"""
# Title: Self-assembly of polycubes (3D and 2D). It stores all complexities of inputs mapping to each output instead of storing the minimum complexity.

# Description: 
This script performs the self-assembly of 2D and 3D polyominoes. 
This is built on top of main_lowKP.py (hydra) but it simplifies the genotype before calculating the complexity.

[Author: Prarthana Agrawal]
[Date (last modified): 10 Aug 2026]
[Version: 0.1]

Older versions:
--------------------------------------------------------------------------------------------------------------------
    # v0.1 (11 feb 2026)
        - Updated complexity calculations (symmetric LZ and minimal stochastic trajectory complexities.)
    # Created using main_lowKP.py: 
        But it simplifies the genotype before calculating the complexity.
    
    # Changes from main.py
        Important distinction from main.py is that this stores all complexities of inputs mapping to each output instead of storing the minimum complexity.
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
import json 
import os


# User-defined packages
#from input import * # input values # use this when input.py is located in the same directory
from core import import_input, get_params, tile_identity_func # tile identities and system parameters
from core.validity_checks import valid_sol_checker # validity of solutions obtained in self-assembly
from utils import get_batch_range, save_all_genotypes, save_assembly_description, get_exh_genotypes
from symmetry import compare_polycubes
#from plots import plot_all_cubes
from utils.genotype_utils import simplify_tiledict, simplify_orientdict
#from utils.binary_utils import tiledict_orientdict_to_binary, lempel_ziv_complexity

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
else: nsplit = 1 #! what does this mean? #! changed this from 0 to 1
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
assembly_settings = {'assembly_type': assembly_type,'Dmax': Dmax,'max_tiles': max_tiles,'kmax': kmax}
parallel_run_config = {'tot_splits': tot_splits,'nsplit': nsplit}

#====================================================================================================================#
#----------------------------------------- Initialisation -----------------------------------------------------------#
#====================================================================================================================#

global UBD, ND, valid
UBD = 0; ND = 0; valid = 0 # number of unbounded (UBD), non deterministic (ND) and valid solutions
sol_stats = dict({'UBD': UBD, 'ND': ND, 'valid': valid}) # solution statistics

#! This is new!
#------------------------- List of dictionaries to store genotypic complexities per phenotype -------------------------#
# Instead of appending complexities values of genotypes mapping to same phenotype to the same list index, we can define dictionaries of each complexity measure.
# Each index will correspond to the shape and the line will be a dictionary with complexity value as keys and number of genotypes mapping to that shape with that complexity value as values. This way we can get the distribution of complexity values for each shape instead of just the minimum complexity value for each shape.
# e.g. old: list = [1,1,1,1,1,1,1
#                 2,2,2,3,4,3,2,2]
#     new: list = [ {1: 7}, 
#                   {2: 5, 3: 2, 4: 1},
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

    #! tile_dict is simplified before complexity calculation for this script. 
    # Simplify tile_dict before assembling
    # 1. '00' all neutral sides 2. '00' sides with no partners 3. renumber sides
    tile_dict = simplify_tiledict(params, tile_dict)
    orient_dict = simplify_orientdict(tile_dict, orient_dict) # simplify orientation dictionary based on simplified tile_dict # assign default orientation wherever tile interface is '00'

    #==============================================================================================================#
    #------------------------------------- Main assembly for this genotype ----------------------------------------#
    #==============================================================================================================#
    
    # check if the resulting phenotype is valid or not
    # NDcomplexity_list (in main_lowKP) was created to return the complexities of all assembly pathways for a genotype, which can then be used to get one complexity per genotype (minimum of the complexities of all assembly pathways for that genotype).
    # Here it is replaced by just complexity (colors)

    output, tile_coord, picked_tiles, complexity, complexity_species, lz_complexity, sol_stats = valid_sol_checker(params, assembly_settings, tile_dict, orient_dict, sol_stats, return_all_complexities=False)

    #! Note: return_all_complexities flag is False because we want one complexity per genotype (minimum of the complexities of all assembly pathways for that genotype). 
    #! Then we can store all complexities for all genotypes mapping to the same shape.
    
    # comp(species) and lz-complexity calculation
    # 1. store complexity along each k-run, 2. take minimum across all k-runs, 3. for a shape, take minimum over all inputs
    #if output == 1: #? Do not compute complexities only for final trajectory, instead compute minima across all k-runs
        #complexity_species = len(set(picked_tiles)) # complexity measure: min number of tile species required

        # Remove tiles with all sides as '00' before calculating lz complexity
        # nonzero_tile_dict = {k: v for k, v in tile_dict.items() if not all(side == '00' for side in v)}
        # nonzero_orient_dict = {k: orient_dict[k] for k in nonzero_tile_dict}
        # binary_genotype = tiledict_orientdict_to_binary(n_sides, nonzero_tile_dict, nonzero_orient_dict)
        # lz_complexity = lempel_ziv_complexity(binary_genotype)
            
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
                    
                    # increment dictionaries (create complexity key if missing)
                    complexity_list[i][complexity] = complexity_list[i].get(complexity, 0) + 1
                    complexity_species_list[i][complexity_species] = complexity_species_list[i].get(complexity_species, 0) + 1
                    lz_complexity_list[i][lz_complexity] = lz_complexity_list[i].get(lz_complexity, 0) + 1

                    #complexity_list[i].append(complexity) # append complexity to the list
                    #complexity_species_list[i].append(complexity_species)
                    #lz_complexity_list[i].append(lz_complexity)
                    shape_type.append(int(shape_type_dict[np.array(valid_shapes[i]).tobytes()]))
                    match = 'yes'
                    break

            if match == 'no':
                valid_shapes.append(tile_coord)
                frequency.append(1)
                complexity_list.append({complexity: 1})
                complexity_species_list.append({complexity_species: 1})
                lz_complexity_list.append({lz_complexity: 1})
                shape_type.append(max(shape_type)+1)
                shape_type_dict[np.array(tile_coord).tobytes()] = max(list(shape_type_dict.values()))+1

        else: # for the first entry
            valid_shapes.append(tile_coord)
            frequency.append(1)
            complexity_list.append({complexity: 1})
            complexity_species_list.append({complexity_species: 1})
            lz_complexity_list.append({lz_complexity: 1})
            shape_type_dict[np.array(tile_coord).tobytes()] = 1
            shape_type.append(int(1))
    else:
        #all_shapes.append(np.array([0])) # 0 implies that the shape is invalid
        shape_type.append(int(0))

plt.show()

# -------------------------------------- Check bookkeeping ---------------------------------------------------------#
# after assembling and before saving
for f_val, hist in zip(frequency, complexity_list):
    assert f_val == sum(hist.values()), "frequency != sum(complexity histogram) — bookkeeping bug"

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
output_dir = "data_files"
os.makedirs(output_dir, exist_ok=True)

file_name = output_dir + "/" + file_name

if dim == 3: save_orientations = True 
else: save_orientations = False

save_all_genotypes(file_name, all_genotypes, all_orientations, save_orientations)
save_assembly_description(file_name, params, assembly_settings, parallel_run_config, master_seed, sol_stats)

np.savetxt(file_name + 'shape_type.txt', shape_type, fmt='%d')
np.savetxt(file_name + 'frequency.txt', frequency, fmt='%d')
#np.savetxt(file_name + 'complexity.txt', complexity_list, fmt='%d')
#np.savetxt(file_name + 'complexity_species.txt', complexity_species_list, fmt='%d')
#np.savetxt(file_name + 'lz_complexity.txt', lz_complexity_list, fmt='%d')

for name, data in zip(['all_complexity', 'all_complexity_species', 'all_lz_complexity'],
                      [complexity_list, complexity_species_list, lz_complexity_list]):
    with open(file_name + name + '.txt', 'w') as f:
        for hist in data:
            # sort by key for deterministic ordering
            sorted_hist = {k: hist[k] for k in sorted(hist)}
            f.write(json.dumps(sorted_hist) + '\n')

if len(valid_shapes) == 1: # if valid_shapes = [[(0,0,0)]], then savetxt interprets it as a 3D array
    with open(file_name + 'valid_shapes.txt', 'w') as f:
        f.write(str(valid_shapes[0]) + '\n')
else:
    valid_shapes_to_save = np.array(valid_shapes, dtype=object)
    np.savetxt(file_name + 'valid_shapes.txt', valid_shapes_to_save, fmt='%s')


# Save minimum, maximum and average complexity for each shape
# Save summary stats (avg, min, max) for each complexity measure
#for name, data in zip(['complexity', 'complexity_species', 'lz_complexity'],
#                      [complexity_list, complexity_species_list, lz_complexity_list]):
#    
    # Compute average, min, max for each shape
#    avg_vals = [sum(d) / len(d) for d in data]
#    max_vals = [max(d) for d in data]
#    min_vals = [min(d) for d in data]
    
#    # Save to files
#    np.savetxt(file_name + name + '_avg.txt', avg_vals, fmt='%.2f')
#    np.savetxt(file_name + name + '_max.txt', max_vals, fmt='%d')
#    np.savetxt(file_name + name + '_min.txt', min_vals, fmt='%d')


