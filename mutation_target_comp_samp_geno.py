""" 
# Title: Group genotypes based on phenotype complexity and mutate them (repeated phenotypes are allowed)

# Objective:
Do genotypes of high complexity phenotypes mutate to genotypes of low complexity phenotypes after mutations? And vice versa? Group the genotypes based on the complexity of their phenotypes and mutate them n_mut times. Analyse how the complexity of the phenotypes change after mutation.
#! Key difference from mutation_target_comp.py: Here we allow repeated phenotypes in each group. In mutation_target_comp.py, we only take unique phenotypes in each group. There, we directly sample from phenotypes. Here, we sample directly from genotypes.

# Description:
    This script analyzes how genotype mutations affect phenotype complexity.
    Steps:
    1. For each genotype, extract its (nsplit, index) and compute three complexity measures:
        - Number of colors
        - Number of species
        - Lempel-Ziv complexity
    2. Store these measures alongside the genotype identifiers in a file.
    3. Sort all genotypes by a chosen complexity metric.
    4. Divide the sorted genotypes into five groups by complexity percentiles:
        - G1: lowest 20%
        - G2: 20-40%
        - G3: 40-60%
        - G4: 60-80%
        - G5: highest 20%
    5. Sample genotypes from each group and apply mutations.
    6. After mutation, record the new complexity measures.
    7. Compare pre- and post-mutation complexities to analyze transitions between complexity groups.
    Note: You can choose to mutate genotypes sequentially (e.g., 2-mut from 1-mut) or independently.

# Notes:
    - You can filter genotypes based on different phenotypic complexity measures:
        1. Number of colors 
        2. Number of species 
        3. Lempel-Ziv complexity 

#-----------------------------------------------------------------------------------------------------#
[Author: Prarthana Agrawal]
[Date: 17 Sep 2025]
[Version: 0.1]
#-----------------------------------------------------------------------------------------------------#
Older versions:
    # built using mutation.py

Example Use:
    python mutation_target_comp_samp_geno.py input.py $complexity_name_for_filter $n_mutations
"""

#=====================================================================================================================#
#--------------------------------------------------- Import modules --------------------------------------------------#
#=====================================================================================================================#
import os
import time
import sys
import random
import pickle
import numpy as np
from tqdm import tqdm # for progress bar

from symmetry import compare_polycubes
from core import import_input, get_params
from core.validity_checks import valid_sol_checker 
from utils import shift_coordinates, save_all_genotypes, save_assembly_description
from utils.binary_utils import  lempel_ziv_complexity, tiledict_orientdict_to_binary
from utils.genotype_utils import extract_genotype, n_point_mutate, get_rep_genotype_indices_for_target_complexity,  convert_genostr_to_tile_dict, simplify_tiledict, get_genotype_groups_based_on_complexity

#=====================================================================================================================#
#----------------------------------------------- Load input.py -------------------------------------------------------#
#=====================================================================================================================#
start = time.time()

# Get the input file path from the command-line argument
if len(sys.argv) < 2:
    print("Please provide the path to input.py")
    sys.exit(1)

input_file_path = sys.argv[1]
input_module = import_input(input_file_path)

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
   
# Parameters
params = get_params(n_tiles, n_sides, neutral_sides, self_int_sides, n_rules, dim)
assembly_settings = {
    'assembly_type': assembly_type,
    'Dmax': Dmax,
    'max_tiles': max_tiles,
    'kmax': kmax
}
parallel_run_config = {'tot_splits': tot_splits}

if dim != 2: raise NotImplementedError('Invalid d. How to mutate orientations? Haven\'t thought of it yet.')

complexity_name_for_filter = 'lz_complexity' # default complexity measure for filtering phenotypes
if len(sys.argv) > 2:
    complexity_name_for_filter = sys.argv[2]

n_mut = 5 # default number of mutations
if len(sys.argv) > 3:
    n_mut = int(sys.argv[3])

#======================================================================================================================#
#------------------------------ Get groups of genotypes and then mutate them ------------------------------------------#
#======================================================================================================================#
path = '' # path to the main directory (not needed when running from subdirectory)
filepath = f'data_files/{n_tiles}s{n_sides}c_combined_' # path to the combined_files

#target_complexity_group = G1, G2, G3, G4, G5 (lowest 20% to highest 20%)
grouped_genotypes_with_complexities = get_genotype_groups_based_on_complexity(path, filepath, tot_splits, complexity_name_for_filter, num_of_groups=5)

for g in range(len(grouped_genotypes_with_complexities)):
    print(f"Number of genotypes in G{g+1}:", len(grouped_genotypes_with_complexities[g]))
    
    grp_genotype_nsplit_index_list = [genotype_tuple[0] for genotype_tuple in grouped_genotypes_with_complexities[g]]
    grp_complexity_list = [genotype_tuple[1] for genotype_tuple in grouped_genotypes_with_complexities[g]]
    grp_complexity_species_list = [genotype_tuple[2] for genotype_tuple in grouped_genotypes_with_complexities[g]]
    grp_lz_complexity_list = [genotype_tuple[3] for genotype_tuple in grouped_genotypes_with_complexities[g]]

    # Save the genotype (nsplit, index) of each group to a file
    mut_file_name = 'mut_files_'+ str(complexity_name_for_filter) + f"_G{g+1}" # folder destination for mutated files
    os.makedirs(mut_file_name, exist_ok=True)
    with open(f"{mut_file_name}/genotype_indices.pkl", "wb") as f:
        pickle.dump(grp_genotype_nsplit_index_list, f)

    #===================================================================================================================#--------------------------------------- Save 0-mut data ----------------------------------------------------#
    #===================================================================================================================

    zero_mut_comp_colors = dict() # 0 mutation: aggregated complexity (colors)
    zero_mut_comp_species = dict() # 0 mutation: aggregated complexity (species)
    zero_mut_lz_comp = dict() # 0 mutation: aggregated lz complexity

    for i in range(len(grp_complexity_list)): # loop over all filtered shapes

        comp_colors = grp_complexity_list[i]
        comp_species = grp_complexity_species_list[i]
        lz_comp = grp_lz_complexity_list[i]

        zero_mut_comp_colors[comp_colors] = zero_mut_comp_colors.get(comp_colors, 0) + 1
        zero_mut_comp_species[comp_species] = zero_mut_comp_species.get(comp_species, 0) + 1
        zero_mut_lz_comp[lz_comp] = zero_mut_lz_comp.get(lz_comp, 0) + 1
        
    filename1 = f"{mut_file_name}/{n_tiles}s{n_sides}c_0mut_agg_" # file name to save the data
    np.savetxt(filename1 + 'complexity.txt', np.array(list(zero_mut_comp_colors.items())), fmt='%s')
    np.savetxt(filename1 + 'complexity_species.txt', np.array(list(zero_mut_comp_species.items())), fmt='%s')
    np.savetxt(filename1 + 'lz_complexity.txt', np.array(list(zero_mut_lz_comp.items())), fmt='%s')

    #===================================================================================================================#--------------------------------- Mutate and save n-mut data ----------------------------------------------------#
    #===================================================================================================================
    for n_mutations in range(1,n_mut+1): # number of mutations to perform

        print("*"*100)
        print(f"Performing {n_mutations} mutations...")

        global UBD, ND, valid
        UBD = 0; ND = 0; valid = 0 # number of unbounded (UBD), non deterministic (ND) and valid solutions
        sol_stats = dict({'UBD': UBD, 'ND': ND, 'valid': valid}) # solution statistics

        # initialise lists
        valid_shapes = list() # list of all "valid" polycubes
        frequency = list() # list of frequency of each valid polycube
        complexity_list = list() # list of complexity of each valid polycube == min number of colors used
        complexity_species_list = list() # alternate complexity measure == min number of species used
        lz_complexity_list = list() # list of Lempel-Ziv complexity of each valid polycube
        shape_type = list() # list of shape "types" (numbered in order of appearance) generated by each genotype
        all_genotypes = list() # list of all original genotypes
        all_mut_genotypes = list() # list of mutated genotypes
        all_orientations = list() # list of all orientations of tiles # not used in this code
        shape_type_dict = dict() # identify a shape to assign a shape "type" to it

        for (nsplit, index) in tqdm(grp_genotype_nsplit_index_list, desc='processing genotypes', unit="genotype", mininterval=300):

            nsplit_filepath = f'data_files/{n_tiles}s{n_sides}c_nsplit_{nsplit}_'
            genotype = extract_genotype(nsplit_filepath, index) # original genotype (returned as str)
            og_tile_dict = convert_genostr_to_tile_dict(genotype) # convert original genotype to tile dictionary
            all_genotypes.append(list(og_tile_dict.values()))
            #print("Original genotype", genotype) 

            #======================================================================================================#
            #------------------ Step 4: Mutate the genotype ---------------------------------------------------# 
            #======================================================================================================#
            # Mutation performed on original 'unsimplified' genotype
            mut_genotype = n_point_mutate(n_sides, dim, genotype, n_mutations, irreducible_mutations=True) #returned as str
            #print("Mutated genotype", mut_genotype)

            #======================================================================================================#
            #--------- Step 5: Convert mutated genotype string to tile_dict and orient_dict ----------------------#
            #======================================================================================================#
            
            tile_dict = convert_genostr_to_tile_dict(mut_genotype) # convert mutated genotype to tile dictionary
            if dim != 2: raise NotImplementedError('only implemented for d=2.')
            orient_dict = {key: 'OOOOUU' for key in tile_dict.keys()}
            all_mut_genotypes.append(list(tile_dict.values()))

            #? New addition
            # Simplify tile_dict before assembling
            # 1. '00' all neutral sides 2. '00' sides with no partners 3. renumber sides
            tile_dict = simplify_tiledict(params, tile_dict)

            #======================================================================================================#
            # ---------- Step 6: Run assembly code for this genotype and get the assembled phenotype. ------------#
            #======================================================================================================#
            
            output, tile_coord, picked_tiles, complexity, sol_stats = valid_sol_checker(params, assembly_settings, tile_dict, orient_dict, sol_stats)

            # comp(species) and lz-complexity calculation
            if output == 1: #? new addition --  no point calculating complexity for invalid shapes
                complexity_species = len(set(picked_tiles)) # complexity measure: min number of tile species required

                # Remove tiles with all sides as '00' before calculating lz complexity
                nonzero_tile_dict = {k: v for k, v in tile_dict.items() if not all(side == '00' for side in v)}
                nonzero_orient_dict = {k: orient_dict[k] for k in nonzero_tile_dict}
                binary_genotype = tiledict_orientdict_to_binary(n_sides, nonzero_tile_dict, nonzero_orient_dict)
                lz_complexity = lempel_ziv_complexity(binary_genotype)
                

            #======================================================================================================#
            #-------------- Step 7: Compare the new phenotype with all old phenotypes ----------------------------#
            #======================================================================================================#

            # default value: the new phenotype is not same as any old encountered phenotype
            match = 'no'
            if output == 1: # shape is valid
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
                shape_type.append(int(0))

        #=============================================================================================================#
        # -------------------- Calculate aggreagated complexity data for each mutation -------------------------------#
        #=============================================================================================================#

        agg_comp_colors_freq = {}
        for comp, freq in zip(complexity_list, frequency):
            agg_comp_colors_freq[comp] = agg_comp_colors_freq.get(comp, 0) + freq

        agg_comp_species_freq = {}
        for comp, freq in zip(complexity_species_list, frequency):
            agg_comp_species_freq[comp] = agg_comp_species_freq.get(comp, 0) + freq

        agg_comp_lz_freq = {}
        for comp, freq in zip(lz_complexity_list, frequency):
            agg_comp_lz_freq[comp] = agg_comp_lz_freq.get(comp, 0) + freq


        filename3 = f"{mut_file_name}/{n_tiles}s{n_sides}c_{n_mutations}mut_agg_" # file name to save the data
        np.savetxt(filename3 + 'complexity.txt', np.array(list(agg_comp_colors_freq.items())), fmt='%s')
        np.savetxt(filename3 + 'complexity_species.txt', np.array(list(agg_comp_species_freq.items())), fmt='%s')
        np.savetxt(filename3 + 'lz_complexity.txt', np.array(list(agg_comp_lz_freq.items())), fmt='%s')

        #==============================================================================================================#
        #------------------- Step 8: Save the data in a file ---------------------------------------------------------#
        #==============================================================================================================#

        file_name = f"{mut_file_name}/{n_tiles}s{n_sides}c_{n_mutations}mut_" # file name to save the data

        if n_mutations == 1: # this is only to store sample genotypes
            save_all_genotypes(f"{mut_file_name}/{n_tiles}s{n_sides}c_0mut_og_", all_genotypes, all_orientations, save_orientations=False)

        save_all_genotypes(file_name + 'mut_', all_mut_genotypes, all_orientations, save_orientations=False)
        master_seed = None # not used
        save_assembly_description(file_name, params, assembly_settings, parallel_run_config, master_seed, sol_stats)

        np.savetxt(file_name + 'shape_type.txt', shape_type, fmt='%d')
        np.savetxt(file_name + 'frequency.txt', frequency, fmt='%d')
        np.savetxt(file_name + 'complexity.txt', complexity_list, fmt='%d')
        np.savetxt(file_name + 'complexity_species.txt', complexity_species_list, fmt='%d')
        np.savetxt(file_name + 'lz_complexity.txt', lz_complexity_list, fmt='%d')

        if len(valid_shapes) == 0:
            valid_shapes_to_save = np.array([],dtype=object) # empty array
            np.savetxt(file_name + 'valid_shapes.txt', valid_shapes_to_save, fmt='%s')

        if len(valid_shapes) == 1: # if valid_shapes = [[(0,0,0)]], then savetxt interprets it as a 3D array
            with open(file_name + 'valid_shapes.txt', 'w') as f:
                f.write(str(valid_shapes[0]) + '\n')
        else:
            valid_shapes_to_save = np.array(valid_shapes, dtype=object)
            np.savetxt(file_name + 'valid_shapes.txt', valid_shapes_to_save, fmt='%s')

print("=" * 100) 
end = time.time()
total_seconds = int(end - start)
hours, remainder = divmod(total_seconds, 3600)
minutes, seconds = divmod(remainder, 60)
print(f"Total run time: {hours}h {minutes}m {seconds}s")