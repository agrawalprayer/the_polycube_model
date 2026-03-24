"""
This combines all parallel runs into a single file and saves it as a .txt file.

Version: 0.6 [17 May 2025]

Older versions:
    # 0.5 [03 Mar 2025]
        - (minor): file destination changed to 'combined_files/' instead of 'data_files/'
        - (minor): updated format of the description file.

    # 0.4 [17 Feb 2025]
        - modification: if lz complexity is empty, don's save it.
        - error: nsplit default value start was set to 0. Changed to 1.
        
    # 0.3 [Date: 14 Feb 2025]
        - modification: added another complexity measure (lempel ziv)
        
    # 0.2 [Date: 28/01/2025]
        - modification: instead of combining all nsplit files into a single file, combine a range of nsplits.
        - shortcoming: very slow to combine all nsplits if tot_splits = 500.

    # 0.1 
        - modification: added another complexity measure (min number of tile species used)

------
# Note: Please run this from the subdirectory where the input.py file is located.
"""

# -------------------------------------- Import packages ----------------------------------------------------------#
import numpy as np
from utils import combine_parallel_runs
from core import import_input, get_params
import sys

#----------------------------- Get input.py path and nsplit range -------------------  -----------------------------#
# Get the input file path from the command-line argument
if len(sys.argv) < 2:
    print("Please provide the path to input.py.")
    sys.exit(1)
input_file_path = sys.argv[1]

#-------------------------------------- Load input.py ----------------------------------------------------------#
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
parallel_run_config = {
    'tot_splits': tot_splits
}

#----------------------------- Get nsplit range -------------------  -----------------------------#
if len(sys.argv) < 4:
    print("Please provide the nsplit ranges. Default is all nsplits.")

if len(sys.argv) > 2:
    start_nsplit = int(sys.argv[2])
    end_nsplit = int(sys.argv[3])
else:
    start_nsplit = 1
    end_nsplit = tot_splits

print('-'*50)
print(f"Combining #{start_nsplit} to #{end_nsplit}.")

#-------------------------------------- Combine parallel runs --------------------------------------------------#
# combine parallel runs
valid_shapes_list, frequency_list, complexity_list, complexity_species_list, lz_complexity_list, sol_stats = combine_parallel_runs(params, parallel_run_config, nsplit_range=(start_nsplit, end_nsplit))
[total_UBD, total_ND, total_valid] = sol_stats

# save combined data files
path = 'combined_files/'

if start_nsplit != 1 or end_nsplit != tot_splits:
    filename = f'{n_tiles}s{n_sides}c_combined_{start_nsplit}to{end_nsplit}nsplits_' # filename prefix
else:
    filename = f'{n_tiles}s{n_sides}c_combined_'

with open(path+filename+'valid_shapes.txt'.format(n_tiles, n_sides), 'w') as file:
    for sublist in valid_shapes_list:
        file.write(f"{sublist}\n")

np.savetxt(path+filename+'frequency.txt'.format(n_tiles, n_sides), frequency_list, fmt='%d')
np.savetxt(path+filename+'complexity.txt'.format(n_tiles, n_sides), complexity_list, fmt='%d')
np.savetxt(path+filename+'complexity_species.txt'.format(n_tiles, n_sides), complexity_species_list, fmt='%d')
if lz_complexity_list:
    np.savetxt(path+filename+'lz_complexity.txt'.format(n_tiles, n_sides), lz_complexity_list, fmt='%d')

# save description file with UBD, ND, valid values
with open(path+filename+'description.txt'.format(n_tiles, n_sides), 'w') as file2:
    file2.writelines(['\n'+'*'*50])
    file2.writelines(['\nNumber of unbounded rules = ', str(total_UBD)])
    file2.writelines(['\nNumber of non deterministic rules = ', str(total_ND)])
    file2.writelines(['\nNumber of valid rules = ', str(total_valid)])
    file2.writelines(['\n'+'*'*50])

pass