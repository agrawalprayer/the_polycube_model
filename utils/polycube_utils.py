"""
# Description: Contains utility functions for polycube self-assembly like return lengths etc

# Functions covered:
    - return_lengths(coord_list)
    - extract_underscore(string, part)
    - shift_coordinates(coord_list)
    - get_bounding_box(coord_list)
    - get_batch_range(params, parallel_run_config)
    - get_exh_genotypes(params)
    - combine_parallel_runs(params, parallel_run_config)
    - get_bar_shift(i, N)
    - zero_singleton_frequencies(frequency_list, shape_content)

# Dependencies:
    - numpy
    - itertools
    - math
    - os

"""

# Import packages
import numpy as np
import itertools
import math

#=============================================================================================================#
def return_lengths(coord_list):
    """ 
    For a passed list of tile coordinates, return lengths and range across x, y and z axes.
    
    Args:
        - coord_list (list of tuples): List of tile origin in the form [(x1, y1, z1), (x2, y2, z2), ...]

    Returns:
        - tuple: lengths of shape across x,y and z axes, minimum values, and maximum values of tile origins 
        (len_x,len_y,len_z), (min_x,min_y,min_z), (max_x,max_y,max_z) 
    """
    coord_array = np.array(coord_list)
    min_vals = np.min(coord_array, axis=0)
    max_vals = np.max(coord_array, axis=0)
    lengths = (max_vals + 1) - min_vals
    return tuple(lengths), tuple(min_vals), tuple(max_vals)

#=============================================================================================================#
def extract_underscore(string, part):
    """ 
    Extracts the first or second part of a string containing digits separated by an underscore '23_4'

    Args:
        - string (str): Input string containing digits separated by an underscore eg '23_4'.
        - part (str): Specifies which part to extract ('first' or 'second').

    Returns:
        - int: Extracted integer value of the specified part.

    Raises:
        - ValueError: If underscore is not found or if the specified part cannot be extracted.
    """
    underscore_index = string.find("_")
    if underscore_index == -1: raise ValueError("Underscore not found in the string.")

    if part == 'first': 
        return int(string[:underscore_index])
    elif part == 'second': 
        if len(string) <= underscore_index+1: 
            raise ValueError("No digit found after underscore.")
        return int(string[underscore_index+1:])
    else:
        raise ValueError("Invalid part specified. Use 'first' or 'second'.")
    
#=============================================================================================================#
def shift_coordinates(coord_list):
    """ 
    Shifts coordinates such that all values are non-negative.

    Args:
        - coord_list (list of tuples): List of coordinates to be shifted in the form [(x1, y1, z1), (x2, y2, z2), ...]

    Returns:
        - shifted_coord_list (list of tuples): Shifted coordinates with all values non-negative.
    """    
    shifted_coord_list = []
    coord_array = np.array(coord_list)
    min_vals = np.min(coord_array, axis=0)
    
    for coord in coord_list:
        shifted_coord = tuple(coord[i] - min_vals[i] for i in range(3))
        shifted_coord_list.append(shifted_coord)

    return shifted_coord_list

#=============================================================================================================#
def get_bounding_box(coord_list):
    """  
    Calculates coordinates of the bounding box of polycube.

    Args:
        - coord_list (list of tuples): list of coordinates of polycube #?(make sure all non-negative values)

    Returns:
        - bbox_coord (list of tuples): list of coordinates of bounding box of polycube

    [Version: 2.0 (25 Feb 2025)]

    #? Big Changes: Changed the bounding box to be a cube with sides equal to the maximum value of x, y, z. 
    #? Earlier bbox could be a cuboid if xmax, ymax, zmax were different. This gave wrong results in comparing
    #? some shapes.
    """
    coord_list = shift_coordinates(coord_list) # make sure all values are non-negative

    coord_array = np.array(coord_list)
    max_vals = np.max(coord_array, axis=0)
    (max_x, max_y, max_z) = max_vals
    
    #! MASSIVE CHANGE 
    #bbox_coord = [
    #    (0,0,0),(max_x+1,0,0),(0,max_y+1,0),(0,0,max_z+1),(max_x+1,max_y+1,0),
    #    (0,max_y+1,max_z+1),(max_x+1,0,max_z+1),(max_x+1,max_y+1,max_z+1)] 
    
    max_val = max(max_x, max_y, max_z)
    bbox_coord = [
        (0,0,0),(max_val+1,0,0),(0,max_val+1,0),(0,0,max_val+1),(max_val+1,max_val+1,0),
        (0,max_val+1,max_val+1),(max_val+1,0,max_val+1),(max_val+1,max_val+1,max_val+1)] 
    # Alternative: not adding +1 to max value creates a bbox surrounding all cube "origins" excluding the outside faces of cubes #! Don't use this because the axes of symmetry don't represent the polycube's axes.
    
    bbox_coord = list(set(bbox_coord))
    return bbox_coord

#==================================================================================================#
def get_batch_range(params, parallel_run_config):
    #! Modify this for more complex parallel runs like not tot_splits = 10
    """
    Calculates the batch range for parallel runs.
    Returns start and end index for each parallel split ("nsplit") in random assembly.
    
    Args:
        - params (dict): {'tile_types':tile_types, 'side_types':side_types, 'neutral_sides':neutral_sides, 'self_int_sides':self_int_sides, 'n_rules':n_rules, 'dim':dim}
        - parallel_run_config (dict): {'tot_splits': tot_splits, 'nsplit': nsplit}

    Returns:
        - indices (list): [start_index, end_index]
            start_index (int): starting index for this nsplit
            end_index (int): last index of genotype range for this nsplit
        - file_name (str): file name to save the data

    eg if total genotypes are 100, nsplits can divide it into batches of 10x10
    start and end index would then be like 0-9, 10-19,..., 90-99
    """
    # load parameter values
    tile_types = params['tile_types']
    side_types = params['side_types']
    n_rules = params['n_rules']
    
    tot_splits = parallel_run_config['tot_splits']
    nsplit = parallel_run_config['nsplit']

    n_tiles = len(tile_types)
    n_sides = len(side_types)

    #----------------------------------- If n_rules is split in parallel runs --------------------------------------#
    if tot_splits > 1:
        split_size = n_rules // tot_splits # usually tot_splits = 10
        start_index = (nsplit - 1) * split_size
        end_index = nsplit * split_size
        if nsplit == tot_splits: # for last split, consider all extra values which could not be divided equally in nsplits
            end_index = n_rules
        file_name = '{}s{}c_nsplit_{}_'.format(n_tiles, n_sides, nsplit)
    else:
        start_index = 0; end_index = n_rules
        file_name = '{}s{}c_'.format(n_tiles, n_sides)

    indices = [start_index, end_index]
    return indices, file_name

#==================================================================================================#

def get_exh_genotypes(params):
    """ 
    [Author: Prarthana Agrawal]
    [Date (last modified): 30 Jan 2025]
    Version: 0.2

    Old Versions:
    # 0.1 
        - modification: automated the dim=2 block using itertools instead of explicit loop.

    Return list of all possible genotypes for exhaustive search.
    Number of genotypes = n_sides**(2*dim*n_tiles)

    Args:
        - params (dict): {'tile_types':tile_types, 'side_types':side_types, 'neutral_sides':neutral_sides, 'self_int_sides':self_int_sides, 'n_rules':n_rules, 'dim':dim}   
            tile_types (list): List of tile types
            side_types (list): List of side types
            neutral_sides (list): List of neutral sides
            self_int_sides (list): List of self-interacting sides
            n_rules (int): Number of rules
            dim (int): Dimension of polycube

    Returns:
        - genotype_combinations_per_tile (list): List of all possible genotypes
    
    Raises:
        - NotImplementedError: If dim is 3.
        - ValueError: If dim is not 2 or 3.

    """

    side_types = params['side_types']
    dim = params['dim']

    if dim == 2: #! test this
        # Generate all possible combinations of 4 sides + 2 neutral placeholders
        genotype_combinations_per_tile = [
            list(side) + ['00', '00'] for side in itertools.product(side_types, repeat=4)
        ]

    #if dim == 2:
    #    genotype_combinations_per_tile = list()
    #    for A1 in side_types:
    #       for A2 in side_types:
    #            for A3 in side_types:
    #                for A4 in side_types:
    #                    side = [A1,A2,A3,A4,'00','00']
    #                    genotype_combinations_per_tile.append(side)

    elif dim == 3:
        raise NotImplementedError("3D exhaustive search not implemented yet.")
        #! Add 3D exhaustive search here if needed
        #! 6 loops over side_types

    else:
        raise ValueError("Invalid dimension. Use 2 or 3.")
    
    return genotype_combinations_per_tile

#==================================================================================================#
def combine_parallel_runs(params, parallel_run_config, nsplit_range):
    """ 
    To combine data from parallel runs into a single file.
    This recalculates the overall frequency, complexity of unique shapes.

    [Author: Prarthana Agrawal]
    [Date (last modified): 04 June 2025]
    Version: 1.0

    Old versions:

    # 0.5 [19 Feb 2025]
        - (major): directly save the first file to lists. No need to compare shapes within the file.
        - (major): for a new file, don't compare with shapes of the same file. Track indices of valid_shapes before adding this file so that redundant comparisons are avoided.
        - (minor): handle single value frequency, complexity, etc. in the file.

    # 0.4 [17 Feb 2025]
        - modification: if lz complexity is not found, skip its calculations.

    # 0.3 [14 Feb 2025]
        - modification: added another complexity measure (lempel ziv)

    # 0.2 [30 Jan 2025]
        - modification: added an option to specify a range of nsplits to combine instead of all nsplits.
        - shortcoming: very slow to combine all nsplits if tot_splits = 500.

    # 0.1 [Jan 2025]
        - modification: added complexity_species to add alternate measure of complexity as number of tile species.
        
    Args:
        - params (dict): {'tile_types':tile_types, 'side_types':side_types, 'neutral_sides':neutral_sides,  'self_int_sides':self_int_sides, 'n_rules':n_rules, 'dim':dim}   
            tile_types (list): List of tile types
            side_types (list): List of side types
            neutral_sides (list): List of neutral sides
            self_int_sides (list): List of self-interacting sides
            n_rules (int): Number of rules
            dim (int): Dimension of polycube
        - parallel_run_config (dict): {'tot_splits': tot_splits, 'nsplit': nsplit}
            tot_splits (int): Total number of parallel splits
            nsplit (int): Number of this split
        - nsplit_range (tuple): Range of nsplits to combine (start_nsplit, end_nsplit)

    Returns:
        - valid_shapes_list (list): List of all valid shape matrices
        - frequency_list (list): List of frequency of each valid shape matrix
        - complexity_list (list): List of complexity of each valid shape
        - complexity_species_list (list): List of complexity (min number of species used)
        - lz_complexity_list (list): List of lempel ziv complexity of each valid shape
        - sol_stats (list): [total_UBD, total_ND, total_valid]
            total_UBD (int): Number of unbounded shapes
            total_ND (int): Number of non-deterministic shapes
            total_valid (int): Number of valid shapes
    
    Raises:
        - FileNotFoundError: If file is missing.    
    """

    import os # for file operations
    from symmetry import compare_polycubes # for comparing shapes
    import re # for regular expressions

    # load parameter values
    tile_types = params['tile_types']
    side_types = params['side_types']
    tot_splits = parallel_run_config['tot_splits']
    n_tiles = len(tile_types)
    n_sides = len(side_types)

    start_nsplit, end_nsplit = nsplit_range

    # ----------------------------------------------------------------------------------------------------------- #
    # Initialise lists with first file values
    # ----------------------------------------------------------------------------------------------------------- #
    first_filepath = f'data_files/{n_tiles}s{n_sides}c_nsplit_{start_nsplit}_'
    frequency_list = np.loadtxt(first_filepath + 'frequency.txt').tolist() # list of frequency of each valid shape matrix
    complexity_list = np.loadtxt(first_filepath + 'complexity.txt').tolist() # list of complexity of each valid shape
    complexity_species_list = np.loadtxt(first_filepath + 'complexity_species.txt').tolist() # list of complexity (min number of species used)
    lz_complexity_list = np.loadtxt(first_filepath + 'lz_complexity.txt').tolist() # list of lempel ziv complexity of each valid shape
    
    # Handle cases where frequency, complexity, etc. are single values
    if np.ndim(frequency_list) == 0: 
        frequency_list = [frequency_list]
        complexity_list = [complexity_list]
        complexity_species_list = [complexity_species_list]
        lz_complexity_list = [lz_complexity_list]

    shape_file = open(first_filepath + 'valid_shapes.txt', 'r')
    shape_content = shape_file.readlines()
    shapes = [eval(line) for line in shape_content]
    valid_shapes_list = list(shapes) # list of all valid shape matrices
    print(valid_shapes_list ) #!delete this

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

    # ----------------------------------------------------------------------------------------------------------- #
    #! deprecated:
    # initialise lists
    #valid_shapes_list = list() # list of all valid shape matrices
    #frequency_list = list() # list of frequency of each valid shape matrix
    #complexity_list = list() # list of complexity of each valid shape
    #complexity_species_list = list() # list of complexity (min number of species used)
    #lz_complexity_list = list() # list of lempel ziv complexity of each valid shape
    #total_UBD = 0 # number of unbounded shapes
    #total_ND = 0 # number of non-deterministic shapes
    #total_valid = 0 # number of valid shapes
    #------------------------------------------------------------------------------------------------------------- #

    lz_complexity_exists = True

    # loop over all parallel runs
    for nsplit in range(start_nsplit+1, end_nsplit+1): # first nsplit is already loaded
        #----------------------------------- Processing split number --------------------------------------#
        print(f"Processing #{nsplit}")
        filepath = 'data_files/{}s{}c_nsplit_{}_'.format(n_tiles, n_sides, nsplit)

        # if file is empty, skip
        if os.path.getsize(filepath + 'frequency.txt') == 0:
            print('empty nsplit file')
            continue

        # load data from each parallel run
        frequency = np.loadtxt(filepath + 'frequency.txt')
        complexity = np.loadtxt(filepath + 'complexity.txt')
        complexity_species = np.loadtxt(filepath + 'complexity_species.txt')
        
        try:
            lz_complexity = np.loadtxt(filepath + 'lz_complexity.txt')
        except OSError:
            lz_complexity_exists = False

        shape_file = open(filepath + 'valid_shapes.txt', 'r')
        shape_content = shape_file.readlines()

        # load shape file
        shapes = []
        for line in shape_content:
            coords = eval(line)
            shapes.append(coords)

        # size of valid shapes file before adding this new file
        max_index_for_comparison = len(valid_shapes_list)

        #! what does this do?
        if frequency.ndim == 0: 
            frequency = [frequency]
            complexity = [complexity]
            complexity_species = [complexity_species]
            if lz_complexity_exists: lz_complexity = [lz_complexity]

            #valid_shapes_list.append(shapes[0])
            #frequency_list.append(float(frequency))
            #complexity_list.append(float(complexity))
            #complexity_species_list.append(float(complexity_species))
            #if lz_complexity_exists:
            #    lz_complexity_list.append(float(lz_complexity))
            #continue

        # Compare shapes and update lists ----------------------------------------------------------------- #
        total_iterations = len(shapes)
        percent_step = 5  # print every 5%
        step_size = max(1, total_iterations * percent_step // 100)  # avoid zero division
        print(f"Total shapes to process: {total_iterations}")

        for i in range(len(shapes)):

            if (i + 1) % step_size == 0 or i == total_iterations - 1:
                percent_complete = ((i + 1) * 100) // total_iterations
                print(f"{percent_complete}% complete ({i + 1}/{total_iterations} shapes processed)")

            match = 'no'
            tile_coord = shapes[i]

            for j in range(max_index_for_comparison): # index of the file being updated
                output1 = compare_polycubes(tile_coord, valid_shapes_list[j]) 
                if output1 == 1: 
                    frequency_list[j] = frequency_list[j] + frequency[i]
                    complexity_list[j] = min(complexity_list[j], complexity[i]) 
                    complexity_species_list[j] = min(complexity_species_list[j], complexity_species[i])
                    if lz_complexity_exists:
                        lz_complexity_list[j] = min(lz_complexity_list[j], lz_complexity[i])
                    match = 'yes'
                    break

            if match == 'no':
                valid_shapes_list.append(tile_coord)
                frequency_list.append(frequency[i])
                complexity_list.append(complexity[i])
                complexity_species_list.append(complexity_species[i])
                if lz_complexity_exists:
                    lz_complexity_list.append(lz_complexity[i])
             
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

    return valid_shapes_list, frequency_list, complexity_list, complexity_species_list, lz_complexity_list, sol_stats

#==================================================================================================#
def get_bar_shift(i, N): 
    """ 
    Shift the bar position for plotting multiple datasets in a bar chart.
    
    Args:
        - i (int): Index of the dataset.
        - N (int): Total number of datasets.

    Returns:
        - float: Shifted position for the bar.
    """
    width = 0.8 / (N + 1)  # bar width/number of datasets
    bar_shift = (i - (N - 1) / 2) * width
    return bar_shift, width

#==================================================================================================#
def zero_singleton_frequencies(frequency_list, shape_content):
    """
    Set frequencies to zero for singleton shape (i.e., shapes with only one cube == [(0,0,0)]).
    This shape corresponds to effectively no assembly.

    Parameters:
        - frequency_list : list[int] or np.ndarray
          Frequencies corresponding to each shape.
        - shape_content : list[str] or list[list[tuple]]
          Given as shape_file.readlines() or list of shapes.
          List of shapes; each shape may be a string or a list of (x, y, z) tuples.

    Returns:
        - new_freq : np.ndarray
            A copy of the frequency list with singleton-shape frequencies set to zero.
    """
    
    new_freq = np.array(frequency_list, dtype=float).copy()

    for i, shape in enumerate(shape_content):
        
        shape_eval = eval(shape.strip())

        # Condition: shape == [(0,0,0)] or length == 1
        if shape_eval == [(0, 0, 0)] or len(shape_eval) == 1:
            new_freq[i] = 0.0

    return new_freq

#==================================================================================================#
def convert_tilecoord_to_2d_matrix(tile_coord):
    """
    Converts tile coordinates of a 2D polycube (with z-coordinate 0) into a 2D matrix representation.
    
    Author: Prarthana Agrawal
    Date: 16 Jan 2026

    Args:
        - tile_coord (list of tuples): List of tile coordinates in the form [(x1, y1, z1), (x2, y2, z2), ...]
    """

    # Check if all z-coordinates in tile_coord are 0
    if not all(z == 0 for _, _, z in tile_coord):
        raise ValueError("All z-coordinates must be 0 for a 2D polycube.")
    
    x_array, y_array, _ = zip(*tile_coord)
    xs, ys = list(x_array), list(y_array)

    min_x, max_x = min(xs), max(xs)
    min_y, max_y = min(ys), max(ys)

    nx = max_x - min_x + 1
    ny = max_y - min_y + 1

    matrix = np.zeros((ny, nx), dtype=int)  # rows = y, cols = x

    for x, y, z in tile_coord:
        row = y - min_y
        col = x - min_x
        matrix[row, col] = 1

    return matrix