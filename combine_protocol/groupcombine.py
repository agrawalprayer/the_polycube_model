"""
This combined files in a heirarchical grouping of 2 till a single combined file is obtained.

[Author: Prarthana Agrawal]
[Version: 0.4 (created on 26 june 2026)]

Difference between combine.py and groupcombine.py is
    - combine.py processes individual nsplit files and combines them (this is level1 of heirarchical grouping)
    - groupcombine.py processes the combined files from combine.py and combines them into a single file (this includes all levels after level1 of heirarchical grouping)

Example: 
    for tot_splits = 10, combine.py will process each nsplit and create level1 combined files:
    1-2 3-4 5-6 7-8 9-10 (after level1 combine.py)
    1-4 5-8 9-10 (combining level1 output using groupcombine.py)
    1-8 9-10 level3
    1-10 level4 is the final output here


Older versions: 
    # 0.3 [Create on 24 June 2026]
        - File2 shapes need be compared in serial order. We can split this into groups and run parallely.
    # 0.2 [Created on: 14 May 2025]
        - Super slow in front loading the first file. Applying cache technique to speed up the first file loading.
    # 0.1 [Date: 17 Feb 2025]
        - earlier known as supercombine.py. This was modified to produce groupcombine.py.

Note: Run it using combine_job.sh script.
"""

# -------------------------------------- Import packages ----------------------------------------------------------#
import numpy as np
import sys
import os
import re
from core import import_input, get_params
from symmetry import compare_polycubes
import psutil
import ast
import time
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed

#-------------------------------------- Load input.py ----------------------------------------------------------#
# Get the input file path from the command-line argument
if len(sys.argv) < 2:
    print("Error: Please provide the path to input.py.")
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

# Chunk size used to divide file2 shapes into groups for parallel processing.
# Defaults to 3000 when parallel_chunk_size is not defined in input.py.
chunk_size = getattr(input_module, "parallel_chunk_size", 3000)

# Parameters
params = get_params(n_tiles, n_sides, neutral_sides, self_int_sides, n_rules, dim)
parallel_run_config = {
    'tot_splits': tot_splits
}

#------------------------------------------------------------------------------------------------------------#
# Track memory usage
def print_memory_usage(note=""):
    process = psutil.Process(os.getpid())
    mem_mb = process.memory_info().rss / 1024 ** 2
    print(f"[MEMORY] {note} - RAM used: {mem_mb:.2f} MB")
#------------------------------------------------------------------------------------------------------------#

#-----------------------------------------------------------------------------------------------------------#
# need to combine files which are already combined in an nsplit range
# for eg combined_1to2.txt, combined_3to4.txt,.. need to combine these files to combined_1to4.txt
#-----------------------------------------------------------------------------------------------------------#

if len(sys.argv) < 3:
    print('Usage: python groupcombine.py input.py "1_2.3_4"')
    sys.exit(1)

range_str = sys.argv[2] # e.g. "1_2.3_4" means combine files: 1to2 and 3to4

try:
    ranges = [tuple(map(int, part.split('_'))) for part in range_str.split('.')]
    if not all(len(r) == 2 for r in ranges):
        raise ValueError("Each range must have exactly two numbers.")
except Exception as e:
    print(f"Error parsing input: {e}")
    sys.exit(1)

#ranges = ((1,160),(161,200)) # manual entry in case of error
print(f"Files to combine: {ranges}", flush=True)
all_values = [value for pair in ranges for value in pair]
min_nsplit = min(all_values)
max_nsplit = max(all_values)

# ----------------------------------------------------------------------------------------------------------- #
# Initialise lists with first file values
# ----------------------------------------------------------------------------------------------------------- #

first_filepath = f'combined_files/{n_tiles}s{n_sides}c_combined_{ranges[0][0]}to{ranges[0][1]}nsplits_'

def load_txt(path, cast=int):
    with open(path, 'r') as f:
        return [cast(line) for line in f if line.strip()]

print("Loading numeric files...", flush=True)
with ThreadPoolExecutor(max_workers=4) as executor:
    f1 = executor.submit(load_txt, first_filepath + 'frequency.txt', int)
    f2 = executor.submit(load_txt, first_filepath + 'complexity.txt', int)
    f3 = executor.submit(load_txt, first_filepath + 'complexity_species.txt', int)
    f4 = executor.submit(load_txt, first_filepath + 'lz_complexity.txt', float)
    frequency_list    = f1.result()
    complexity_list   = f2.result()
    complexity_species_list = f3.result()
    lz_complexity_list = f4.result()

print("Loading shapes file...", flush=True)

t = time.time()
# with open(first_filepath + 'valid_shapes.txt', 'r') as shape_file:
#     raw_lines = [line for line in shape_file if line.strip()]
# print("read time:", time.time() - t)

# t = time.time()
# valid_shapes_list = [ast.literal_eval(line) for line in raw_lines]

with open(first_filepath + 'valid_shapes.txt', 'r') as shape_file:
    valid_shapes_list = [ast.literal_eval(line) for line in shape_file if line.strip()]
print("parse time:", time.time() - t)

print_memory_usage("After loading first file")

with open(first_filepath +'description.txt', 'r') as file1:
    content = file1.read()

    # Use regular expressions to find the relevant numbers
    ubd_match = re.search(r'Number of unbounded rules = (\d+)', content)
    nd_match = re.search(r'Number of non deterministic rules = (\d+)', content)
    valid_match = re.search(r'Number of valid rules = (\d+)', content)
    
    # Extract values and return them as integers
    total_UBD = int(ubd_match.group(1)) if ubd_match else 0
    total_ND = int(nd_match.group(1)) if nd_match else 0
    total_valid = int(valid_match.group(1)) if valid_match else 0
    
#-------------------------------------------------------------------------------------------------------#
# NEW: worker function for one chunk of file-2 shapes
def process_shape_group(args):
    group_id, idx_list, shapes, frequency, complexity, complexity_species, lz_complexity, ref_shapes = args

    # NEW: local updates to existing file-1 shapes
    matched_updates = {}

    # NEW: shapes in file 2 that do not match any file-1 shape
    new_shapes = []

    for i in idx_list:
        tile_coord = shapes[i]
        matched_index = -1

        for j, ref_shape in enumerate(ref_shapes):
            if compare_polycubes(tile_coord, ref_shape) == 1:
                matched_index = j
                break

        if matched_index >= 0:
            if matched_index not in matched_updates:
                matched_updates[matched_index] = [
                    int(frequency[i]),
                    int(complexity[i]),
                    int(complexity_species[i]),
                    float(lz_complexity[i]),
                ]
            else:
                matched_updates[matched_index][0] += int(frequency[i])
                matched_updates[matched_index][1] = min(matched_updates[matched_index][1], int(complexity[i]))
                matched_updates[matched_index][2] = min(matched_updates[matched_index][2], int(complexity_species[i]))
                matched_updates[matched_index][3] = min(matched_updates[matched_index][3], float(lz_complexity[i]))
        else:
            new_shapes.append((
                tile_coord,
                int(frequency[i]),
                int(complexity[i]),
                int(complexity_species[i]),
                float(lz_complexity[i]),
            ))

    return group_id, matched_updates, new_shapes

# loop over all parallel runs
for (start_nsplit, end_nsplit) in ranges[1:]: # from second pair onwards
    #----------------------------------- Processing split number --------------------------------------#
    print(f"Processing split range: #{start_nsplit}-#{end_nsplit}")
    filepath = 'combined_files/{}s{}c_combined_{}to{}nsplits_'.format(n_tiles, n_sides, start_nsplit, end_nsplit)

    # check if the file exists
    if not os.path.exists(filepath + 'frequency.txt'):
        print(f"File {filepath + 'frequency.txt'} does not exist. Stopping execution.")
        sys.exit(1)

    # if file is empty, skip
    if os.path.getsize(filepath + 'frequency.txt') == 0:
        print('empty nsplit file')
        continue

    # load data from each parallel run
    with ThreadPoolExecutor(max_workers=4) as executor:
        f1 = executor.submit(load_txt, filepath + 'frequency.txt', int)
        f2 = executor.submit(load_txt, filepath + 'complexity.txt', int)
        f3 = executor.submit(load_txt, filepath + 'complexity_species.txt', int)
        f4 = executor.submit(load_txt, filepath + 'lz_complexity.txt', float)
        frequency    = np.array(f1.result())
        complexity   = np.array(f2.result())
        complexity_species = np.array(f3.result())
        lz_complexity = np.array(f4.result())

    with open(filepath + 'valid_shapes.txt', 'r') as shape_file:
        shapes = [ast.literal_eval(line) for line in shape_file if line.strip()]

    print_memory_usage("After loading new file")

    # size of valid shapes file before adding this new file
    # UPDATED: size of file-1 shapes before adding file-2 shapes
    max_index_for_comparison = len(valid_shapes_list)

    # NEW: choose worker count from Slurm if available
    slurm_cpus = int(os.environ.get("SLURM_CPUS_PER_TASK", os.cpu_count() or 1))

    # NEW: split file-2 shapes into chunks of 3000 shapes each
    #chunk_size = 1000
    n_groups = max(1, (len(shapes) + chunk_size - 1) // chunk_size)

    # NEW: do not use more workers than groups
    n_workers = min(slurm_cpus, n_groups)

    print(f"Total shapes to process: {len(shapes)}", flush=True)
    print(f"Using {n_groups} groups and {n_workers} workers", flush=True)

    # NEW: split by indices, not by copying shapes
    shape_index_groups = np.array_split(np.arange(len(shapes)), n_groups)

    tasks = [
        (
            group_id,
            idx_group.tolist(),
            shapes,
            frequency,
            complexity,
            complexity_species,
            lz_complexity,
            valid_shapes_list,   # file-1 shapes are read-only here
        )
        for group_id, idx_group in enumerate(shape_index_groups)
        if len(idx_group) > 0
    ]

    # NEW: run chunks in parallel
    group_results = []
    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = {
            executor.submit(process_shape_group, task): task[0]
            for task in tasks
        }

        for fut in as_completed(futures):
            group_id = futures[fut]
            print(f"Group {group_id} completed.", flush=True)

            group_results.append(fut.result())

    # NEW: merge deterministically after all workers finish
    group_results.sort(key=lambda x: x[0])

    for group_id, matched_updates, new_shapes in group_results:
        # NEW: apply frequency/complexity updates to existing file-1 shapes
        for j, upd in matched_updates.items():
            freq_add, comp_new, comp_species_new, lz_new = upd
            frequency_list[j] += freq_add
            complexity_list[j] = min(complexity_list[j], comp_new)
            complexity_species_list[j] = min(complexity_species_list[j], comp_species_new)
            lz_complexity_list[j] = min(lz_complexity_list[j], lz_new)

        # NEW: append shapes that were not found in file 1
        for tile_coord, freq_val, comp_val, comp_species_val, lz_val in new_shapes:
            valid_shapes_list.append(tile_coord)
            frequency_list.append(freq_val)
            complexity_list.append(comp_val)
            complexity_species_list.append(comp_species_val)
            lz_complexity_list.append(lz_val)

            
    # Combine description files
    with open(filepath +'description.txt', 'r') as file1:
        content = file1.read()
    
        # Use regular expressions to find the relevant numbers
        ubd_match = re.search(r'Number of unbounded rules = (\d+)', content)
        nd_match = re.search(r'Number of non deterministic rules = (\d+)', content)
        valid_match = re.search(r'Number of valid rules = (\d+)', content)
        
        # Extract values and return them as integers
        ubd = int(ubd_match.group(1)) if ubd_match else 0
        nd = int(nd_match.group(1)) if nd_match else 0
        valid = int(valid_match.group(1)) if valid_match else 0
        
    # Update total values
    total_UBD += ubd
    total_ND += nd
    total_valid += valid

    sol_stats = [total_UBD, total_ND, total_valid]

#-------------------------------------- Save combined data files --------------------------------------------------#
path = 'combined_files/'
filename = f'{n_tiles}s{n_sides}c_combined_{min_nsplit}to{max_nsplit}nsplits_'

with open(path+filename+'valid_shapes.txt'.format(n_tiles, n_sides), 'w') as file:
    for sublist in valid_shapes_list:
        file.write(f"{sublist}\n")

np.savetxt(path+filename+'frequency.txt'.format(n_tiles, n_sides), frequency_list, fmt='%d')
np.savetxt(path+filename+'complexity.txt'.format(n_tiles, n_sides), complexity_list, fmt='%d')
np.savetxt(path+filename+'complexity_species.txt'.format(n_tiles, n_sides), complexity_species_list, fmt='%d')
np.savetxt(path+filename+'lz_complexity.txt'.format(n_tiles, n_sides), lz_complexity_list, fmt='%.3f')

# save description file with UBD, ND, valid values
with open(path+filename+'description.txt'.format(n_tiles, n_sides), 'w') as file2:
    file2.writelines(['\n'+'*'*50])
    file2.writelines(['\nNumber of unbounded rules = ', str(total_UBD)])
    file2.writelines(['\nNumber of non deterministic rules = ', str(total_ND)])
    file2.writelines(['\nNumber of valid rules = ', str(total_valid)])
    file2.writelines(['\n'+'*'*50])

pass

print_memory_usage("After completion:")