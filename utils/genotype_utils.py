""" 
Title: Utility functions to extract genotypes, mutate them or perform other calculations.
[Date: 29 Jan 2025]
[Author: Prarthana Agrawal]

Description:
    This file contains functions to perform various analysis on genotypes.

Functions:
    - get_sample_genotype_indices(path, filepath, n_tiles, n_sides, tot_splits, sample_size): to find representative genotype positions in the form of (nsplit, index) for each valid phenotype
    - pick_random_index_for_shape(args): to randomly pick an index for a given shape in the nsplit_shape_type file.
    - extract_genotype(full_path, index): get genotype at a specific nsplit (indicated in full_path) and index
    - extract_orientation(full_path, index): get orientations at a specific nsplit (indicated in full_path) and index
    - n_point_mutate(n_sides, dim, genotype, n_mutations=1, irreducible_mutations=False): to point mutate n_mutations times a genotype and return mutated genotype.
    - locate_genotype(path, filepath, n_tiles, n_sides, genotype, tot_splits=10): locate genotype in the all_genotypes file.
    - simplify_genotype(params, flat_genotype, flat_orientation): simplify the genotype by zeroing neutral sides, sides with no partners, and removing tiles with all sides '00'.
    - convert_genostr_to_tile_dict(genotype_string): convert string of genotype (list of lists written in a file as string) to tile_dict.
    - convert_genostr_to_genolist(genostr): convert genotype string (read from the file) to list of lists.
    - flatten_dict(tile_dict): flatten the tile_dict to a flat list of sides.
    - unflatten_list(flat_list): convert a flat list of sides back into a dictionary format with tile names A-Z.
    - simplify_tiledict(params, tile_dict): simplify the given tile_dict under three simplification rules
    - get_genotype_for_this_shape(args): get the genotype for a specific shape.
    - get_rep_genotype_indices_for_target_complexity(args): to find representative genotype positions in the form of (nsplit, index) for each valid phenotype belonging to each group of complexity (G1 to G5).

#Note: This file was previously located in analyzer. Moved to utils (3 mar 2025)

Dependencies:
    - random
    - linecache
    - string
    - utils.shift_coordinates
    - symmetry.compare_polycubes
    - plots.plot_all_cubes
    - core.simplify_genotype
    - utils.rqrd_num_bits
    - utils.number_to_binary

"""

# Load packages
import numpy as np
import random
import linecache
import string
from more_itertools import sort_together

from utils import shift_coordinates
from symmetry import compare_polycubes
from plots import plot_all_cubes

def get_sample_genotype_indices(path, filepath, n_tiles, n_sides, tot_splits, sample_size, num_of_cores=1):
    """ 
    Get representative genotype for each unique phenotype.

    [Author: Prarthana Agrawal]
    [Version: 0.3 (21 May 2025)]

    # Older versions:

        - 0.2 (20 May 2025)
        (major): parallelized the process of finding (nsplit, index) values for each valid shape in the combined file.

        - 0.1 
        (major): instead of first storing all (nsplit, index) values and then randomly sampling from them,
        we now randomly sample the indices while iterating through the nsplit_shape_type file picked at random.

    Args:
        - path (str): path to the folder where data files live.
        - filepath (str): folder/file details within the main directory.
        - n_tiles (int): Number of tile species.
        - n_sides (int): Number of sides/ colors.
        - tot_splits (int): Number of splits for the parallel run.
        - sample_size (int): Number of sample genotypes to be generated for each phenotype.
    
    Returns:
        - sample_indices_dict (dict): Dictionary with key as the original shape number in combined file and value as the list of sample indices (tuple of nsplit and index for each phenotype.
        eg {1: [(2,1032), (3,43)], 2: [(8,31)]}
    """

    import multiprocessing

    # Step1: Load combined valid shapes data file and extract the final valid shapes.
    # Step2: For each unique shape, find its corresponding number in each nsplit_valid_shapes.
    # Step3: Randomly pick nsplits from the list of nsplits where the shape is present.
    # Step4: For each nsplit, pick an index where the shape is found.
    # Step5: Return the list of (nsplit, index) values for each shape.

    # ==========================================================================================================# 
    # ---------------------------------- Load combined valid shapes data -------------------------------------- #
    # ==========================================================================================================# 
    """
    Step1: Load combined_valid_shapes data file and extract the final valid shapes.
    """
    shape_file = open(path + filepath + 'valid_shapes.txt', 'r')
    shape_content = shape_file.readlines()

    combined_shapes = [] # final valid shapes from combined file
    for line in shape_content:
        coords = eval(line.strip())
        shift_coords = shift_coordinates(coords) # make coordinates positive
        combined_shapes.append(shift_coords)

    #print("Unique shapes \n", combined_shapes)
    #print("Number of unique shapes: ", len(combined_shapes))

    # ==========================================================================================================# 
    # ------------For each unique shape, find its corresponding number in each nsplit_valid_shapes ------------ #
    # ==========================================================================================================# 
    """
    Step2: For each unique shape, find its corresponding number in each nsplit_valid_shapes.
    For example [(0,0,0)] can be shape number #1 in combined_valid_shapes but #3 in nsplit=2.
    This is because in each nsplit, the shapes are numbered in order of appearance.
    """

    shape_num_dict = dict() # key: shape number, value: list of corresponding shape numbers in each nsplit
    sample_indices_dict = dict() # key: original shape number in combined file, value: list of sample indices

    for num1, shape1 in enumerate(combined_shapes,1):
        #print("**************** Shape number: ", num1, "****************")
        match_list = list() # list of corresponding shape numbers in each nsplit

        for nsplit in range(1, tot_splits+1):
            
            filepath = f"data_files/{n_tiles}s{n_sides}c_nsplit_{nsplit}_"
            nsplit_shape_file = open(path + filepath + 'valid_shapes.txt', 'r')
            nsplit_shape_content = nsplit_shape_file.readlines()
            nsplit_shapes = [eval(line.strip()) for line in nsplit_shape_content] # list of shapes in nsplit

            found = 'no'
            for num2, shape2 in enumerate(nsplit_shapes,1):
                if compare_polycubes(shape1, shape2) == 1: #! make sure the function shifts coordinates first
                    """
                    %matplotlib inline
                    params = dict()
                    fig = plt.figure(figsize=(8, 6))
                    
                    ax1 = fig.add_subplot(121, projection='3d')
                    ax1.view_init(elev=-149, azim=138)
                    plot_all_cubes(params, ax1, shape1, picked_tiles=[], cube_outline='False', axes_lines='False')
                    ax1.set_title(f"Shape {num1}")

                    ax2 = fig.add_subplot(122, projection='3d')
                    ax2.view_init(elev=-149, azim=138)
                    plot_all_cubes(params, ax2, shape2, picked_tiles=[], cube_outline='False', axes_lines='False')
                    ax2.set_title(f"Shape {num2}")
                    
                    plt.show()
                    """
                    found = 'yes'
                    match_list.append(num2) # shape 'type' number in nsplit
                    break

            if found == 'no':
                match_list.append(0)

        #print(f"Shape {num1} is represented by: \n", match_list)
        shape_num_dict[num1] = match_list
        #print("Match list: ", match_list)
        
        #----------------------------------------------------------------------------------------------------------#
        # Filter nsplits where the shape is actually present, along with its shape number
        nsplit_matchnum_list = list() # list of tuples (nsplit, shape number in nsplit)
        valid_nsplits = [split for split in range(1, tot_splits+1) if match_list[split-1] != 0]

        for split in valid_nsplits:
            match_num = match_list[split-1]
            nsplit_matchnum_list.append((split, match_num))
        
        #print("(nsplit, shape number in nsplit) \n", nsplit_matchnum_list)

        #==========================================================================================================#    #-----------Search shape_type file and randomly pick indices where desired shape type exists---------------#
        #==========================================================================================================#
        """
        Step3: Randomly pick nsplits from the list of nsplits where the shape is present.
        For each nsplit, randomly pick an index until the shape is found.
        Repeat for sample_size times.

        Step3 (deprecated): Open nsplit_shape_type files and find (nsplit, index) values where the shape type if found.
        we now have a list of nsplits where the shape is present and the corresponding shape number for that nsplit.
        we can then later extract genotypes at this (nsplit, index) values.
        """

        # randomly pick nsplit from the list of nsplits where the shape is present
        random_nsplit_w_shapenum_list = [random.choice(nsplit_matchnum_list) for _ in range(sample_size)]

        # Parallel reservoir sampling
        with multiprocessing.Pool(processes=num_of_cores) as pool:
            args_list = [(path, n_tiles, n_sides, nsplit, match_num) for (nsplit, match_num) in random_nsplit_w_shapenum_list]
            # parallelize the process of finding (nsplit, index) values for each valid shape in the combined file
            sample_indices = pool.map(pick_random_index_for_shape, args_list) 

        sample_indices_dict[num1] = sample_indices

    return sample_indices_dict


#==========================================================================================================#
def pick_random_index_for_shape(args):
    """ 
    For the passed nsplit, find indices where the searched match_num shape is present and then randomly pick one of those indices. This employs reservoir sampling.
    This is Step 4 of the full process.
    This is done in parallel for all shapes.

    Args:
        - args (tuple): Tuple containing path, n_tiles, n_sides, nsplit, match_num.
    
    Returns:
        - (nsplit, chosen_index) (tuple): Tuple containing nsplit and the randomly chosen matching index.
    """
    path, n_tiles, n_sides, nsplit, match_num = args
    nsplit_filepath = f"data_files/{n_tiles}s{n_sides}c_nsplit_{nsplit}_" # path to the nsplit data files
    shape_type_path = path + nsplit_filepath + 'shape_type.txt'

    # Open the shape_type file for the given nsplit.
    # Find indices where the shape type is equal to match_num.
    # Use reservoir sampling to randomly pick one of those indices.

    with open(shape_type_path, 'r') as f:
        count = 0 
        chosen_index = None

        for idx, shape_type in enumerate(f):
            if int(shape_type) == match_num:
                count += 1
                if random.randint(1, count) == 1:
                    chosen_index = idx

    if chosen_index is None:
        raise ValueError(f"Shape number {match_num} not found in nsplit {nsplit} shape_type file.")

    # For the passed nsplit and shape_num, return the randomly chosen index in the nsplit file where the shape is present. This can help to choose genotype at this index.
    return (nsplit, chosen_index)

#==========================================================================================================# 
def extract_genotype(full_path, index):
    """
    Extract genotype from the data file at given nsplit and index. 
    [Author: Prarthana Agrawal]
    [Date: 27 Jan 2025]

    Args:
        - full_path (str): full path of the genotypes file
            # Note: full_path should include correct nsplit details
        - index (int): Index at which the genotype is present in the data file.
        Note: line number = index + 1

    Returns:
        - genotype (str): string of genotype (as list of lists)
    """
    
    line_num = index+1
    line = linecache.getline(full_path + 'all_genotypes.txt', line_num)
    if line: genotype = line.strip()
    else: raise ValueError(f"Line {line_num} does not exist in the file.")
    
    #genotype_as_list = [lst.strip("[]").replace("'", "").split(", ") for lst in genotype.split("] [")]

    return genotype


#=====================================================================================================#
def extract_orientation(full_path, index):
    """
    Extract orientations from the data file at given nsplit and index. 
    [Author: Prarthana Agrawal]
    [Date: 12 Feb 2025]

    Args:
        - full_path (str): full path of the orientation file
            # Note: full_path should include correct nsplit details
        - index (int): Index at which the orientation is present in the data file.
        Note: line number = index + 1

    Returns:
        - orientations (list): orientations as list of lists
    """

    line_num = index+1
    line = linecache.getline(full_path + 'all_orientations.txt', line_num)
    if line: orientations = line.strip()
    else: raise ValueError(f"Line {line_num} does not exist in the file.")
    
    orientations_as_list = orientations.split()

    return orientations_as_list

#====================================================================================================#

def n_point_mutate(n_sides, dim, tile_dict, orient_dict=None, n_mutations=1, irreducible_mutations=False,       mutate_what='both'):
    """ 
    Point mutate a genotype and return the mutated genotype.
    For 2D, only interface types are mutated.
    For 3D, both interface types and orientations can be mutated.

    Author: Prarthana Agrawal
    Version: 1.0 (19 Jan 2026)

    Older version:
        - major change over v0.1
          extended to 3D polycubes where orientations can also be mutated.
          input genotype is as tile_dict now instead of a string. 
          #! If genotype is input from the txt file, first convert to a tile dict format and then input to this function

    Args:
        - n_sides (int): Number of sides/ colors (to know what to mutate into).
        - dim (int): dimensionality of the assembly (2 or 3).
        - tile_dict (dict): Original tile interfaces in tile_dict dictionary format.
        - orient_dict (dict): Original tile orientations in orient_dict dictionary format.
        - n_mutations (int): Number of point mutations to perform. (default value = 1)
        - irreducible_mutations (bool): True means any nth mutation is not reducible to (n-1)th mutation.
            eg 2-mut should not be an effective 1-mut meaning that the same interface was mutated twice.
        - mutate_what (str): 'interface', 'orientation', or 'both' to specify what to mutate.

        (deprecated)genotype (str): Original genotype read straight from the all_genotypes file. (This is the tile_dict)

    Returns:
        - mut_tile_dict (dict): Mutated tile interfaces in tile_dict format.
        - mut_orient_dict (dict): Mutated orientations in orient_dict format.
        (deprecated) mut_genotype (str): Mutated genotype in same string format as original file.

    #! Note: here 2-mutations may or may not mean that it can be an effective 1-mut. Choose what you want!
    """    
    # ==================================================================================================#
    # Input validations
    # ==================================================================================================#
    if dim == 3 and orient_dict is None:
        raise ValueError("For 3D mutations, orient_dict must be provided")
    
    if dim not in [2,3]:
        raise ValueError("Dimension must be either 2 or 3")
    
    # ==================================================================================================#
    # Orientation and interface alternatives for mutation
    # ==================================================================================================#
    # Interfaces
    side_types = list(range(0,n_sides)) # types of sides/interfaces 
    formatter = "{:02d}".format
    side_types = list(map(formatter,side_types))
    
    # Orientations
    def orient_types_available(side_index):
        if side_index in [0,2]: orient_types = 'IOLR' # North-South options: In, out, left, right
        elif side_index in [1,3]: orient_types = 'IOUD' # East-West options: In, out, up, down
        elif side_index in [4,5]: orient_types = 'UDLR' # Back-Front options: up, down, left, right
        else: raise ValueError("Invalid side index")
        return orient_types

    # ==================================================================================================#
    # Convert tile and orient dicts to list format for easier manipulation
    # ==================================================================================================#
    interface_list = [tile_dict[tile] for tile in tile_dict.keys()]
    if dim == 3: orient_list = [orient_dict[tile] for tile in orient_dict.keys()]

    # Truncate the last two interfaces (they are '00' by default) for 2D
    if dim == 2:
        if any(lst[-2:] != ['00', '00'] for lst in interface_list): 
            raise ValueError("For 2D, the last two interfaces should be zero.")
        interface_list = [lst[:-2] for lst in interface_list]

    # Copy original lists to mutated lists
    mut_interface_list = [lst.copy() for lst in interface_list]
    if dim == 3: 
        mut_orient_list = [orient_str for orient_str in orient_list]
        if len(mut_interface_list) != len(mut_orient_list): 
            raise ValueError("Interface list and orientation list lengths do not match")
        
    # ==================================================================================================#
    # Perform mutations
    # ==================================================================================================#
    mutations_done = 0
    mutated_positions = set()  # to keep track of mutated positions and choice of interface/orientation
    # stores positions as (tile index, side position, 'interface'/'orientation')

    while mutations_done < n_mutations:

        # Pick a tile index and side position
        # here, using mut_interface_list does not matter since both lists have same length
        tile_idx = random.randint(0, len(mut_interface_list) - 1)  # choose a sublist to mutate -- which tile to mutate
        interface_sublist = mut_interface_list[tile_idx]  # choose a sublist to mutate -- at the random index
        pos = random.randint(0, len(interface_sublist) - 1)  # choose a position to mutate -- which side to mutate
        
        # Determine whether to mutate interface or orientation
        if dim == 2: interface_or_orient = 'interface'  # only interface can be mutated in 2D
        elif dim == 3: 
            if mutate_what == 'interface': interface_or_orient = 'interface'
            elif mutate_what == 'orientation': interface_or_orient = 'orientation'
            elif mutate_what == 'both':
                interface_or_orient = random.choice(['interface', 'orientation'])
            else: raise ValueError("Invalid mutate_what value")

        # Check if irreducible_mutations is enabled and the position has been mutated before
        if irreducible_mutations and (tile_idx, pos, interface_or_orient) in mutated_positions:
            continue  # skip this mutation and try again

        # Mutate interface if chosen
        if interface_or_orient == 'interface':
            original_val = interface_sublist[pos]  # original interface at the position
            new_val = random.choice([val for val in side_types if val != original_val])  # new mutated interface
            interface_sublist[pos] = new_val  # update the tile with the new interface
            # here interface_sublist is alias for mut_interface_list[tile_idx], so no need to reassign

        # Mutate orientation if chosen
        elif interface_or_orient == 'orientation':
            orient_string = mut_orient_list[tile_idx]  # get the orientation sublist as a list

            if len(orient_string) != len(interface_sublist):
                raise ValueError("Orientation string length does not match interface sublist length.")
            
            original_orient = orient_string[pos]  # original orientation at the position
            
            if original_orient not in orient_types_available(side_index=pos):
                raise ValueError("Something's wrong. Original orientation not permissible for this side position=", pos)
            
            available_orients = orient_types_available(side_index=pos)
            new_orient = random.choice([o for o in available_orients if o != original_orient])# new mutated orientation
            orient_string = orient_string[:pos] + new_orient + orient_string[pos+1:]  # update the orientation string
            mut_orient_list[tile_idx] = orient_string  # update the mutated orientation list

        else: raise ValueError("Invalid choice for interface_or_orient")

        mutations_done += 1
        mutated_positions.add((tile_idx, pos, interface_or_orient))  # add the position to the set of mutated positions
             
    if dim == 2:
        # if 2D and lists were truncated, then add back the '00' interfaces to be consistent with 3D calculations
        [lst.extend(['00', '00']) for lst in mut_interface_list]
    
    #====================================================================================================#
    # Convert mutated interface and orientation lists back to dict format
    #====================================================================================================#
    mut_tile_dict = {}
    mut_orient_dict = {}
    tile_keys = list(tile_dict.keys())

    for i, tile_name in enumerate(tile_keys):
        mut_tile_dict[tile_name] = mut_interface_list[i]
        if dim == 3:
            mut_orient_dict[tile_name] = mut_orient_list[i]
    
    if dim == 2: mut_orient_dict = orient_dict  # for 2D, return the original orient_dict as is
    return mut_tile_dict, mut_orient_dict
    
#====================================================================================================#
def locate_genotype(path, filepath, n_tiles, n_sides, genotype, tot_splits=10):
    """ 
    Given a genotype string, locate its position in the all_genotypes file.

    Args:
        - path (str): path to the folder where data files live.
        - filepath (str): folder/file details within the main directory
        - n_tiles (int): Number of tile species.
        - n_sides (int): Number of sides/ colors.
        - genotype (str): Genotype string to locate.
        - tot_splits (int): Number of splits for the parallel run.

    Returns:
        - (nsplit, line number) (tuple): Position of the genotype in the all_genotypes file.
        Note: line number = index + 1
    
    """

    for nsplit in range(1, tot_splits+1):
        filepath = f"data_files/{n_tiles}s{n_sides}c_nsplit_{nsplit}_"
    
        with open(path + filepath + 'all_genotypes.txt', 'r') as file:
            for linenum, line in enumerate(file, start=1):
                if genotype in line:
                    #print(f"Genotype found at line number {linenum} in nsplit {nsplit}")
                    return (nsplit, linenum)
                    
#====================================================================================================#
def simplify_genotype(params, flat_genotype, flat_orientation):
    """ 
    Simplifies the genotype by zeroing neutral sides, sides with no partners, and removing tiles with all sides '00'.
    Renumber interfaces respecting their rules of attachment.

    [Author: Prarthana Agrawal]
    [Date: 13 Feb 2025]

    Args:
        - params (dict): {'tile_types':tile_types, 'side_types':side_types, 'neutral_sides':neutral_sides, 'self_int_sides':self_int_sides, 'n_rules':n_rules, 'dim':dim}
        - flat_genotype (list): List of all sides of all tiles.
        - flat_orientation (list): List of all orientations of all tiles.

    Returns:
        - flat_genotype (list): Simplified genotype.
        - flat_orientation (list): Simplified corresponding orientations.
    """
    
    from core import rules_func

    # Extract parameters
    neutral_sides = params['neutral_sides']
    
    # Load rules of attachment
    rules_dict = rules_func(params)

    #print("Given genotype: ", flat_genotype)
    #print("Given orientation: ", flat_orientation)

    # Simplification #1 ---------------------------------------------------------------
    # Zero all neutral sides
    flat_genotype = ['00' if side in neutral_sides else side for side in flat_genotype]
    #print("All sides after zeroing all neutral sides: ", flat_genotype)

    # Simplification #2 ---------------------------------------------------------------
    # Zero all sides with no partners
    flat_genotype = [side if (side == '00' or rules_dict[side] in flat_genotype) else '00'
    for side in flat_genotype]
    #print("All sides after zeroing all sides with no partners: ", flat_genotype)

    # Simplification #3 ---------------------------------------------------------------
    # If a tile_dict key has all values '00', remove that key
    # Splice flat_genotype into sublists of 6 consecutive elements each
    flat_genotype_sublists = [flat_genotype[i:i + 6] for i in range(0, len(flat_genotype), 6)]
    
    if flat_orientation is not None:
        flat_orientation_sublists = [flat_orientation[i:i + 6] for i in range(0, len(flat_orientation), 6)]
        #print("Flat genotype sublists: ", flat_genotype_sublists)
        #print("Flat orientation sublists: ", flat_orientation_sublists)

        # Remove sublists where all 6 elements are '00' and corresponding orientation sublists
        filtered_genotype_orientation = [(genotype_sublist, orientation_sublist) for genotype_sublist, orientation_sublist in zip(flat_genotype_sublists, flat_orientation_sublists) if not all(side == '00' for side in genotype_sublist)]
        #print("Filtered genotype and orientation: ", filtered_genotype_orientation)
        
        flat_genotype_sublists, flat_orientation_sublists = zip(*filtered_genotype_orientation) if filtered_genotype_orientation else ([], [])
        #print("Flat genotype sublists after removing all '00' sublists: ", flat_genotype_sublists)

        # Flatten the list back to a single list
        flat_genotype = [side for sublist in flat_genotype_sublists for side in sublist]
        flat_orientation = [side for sublist in flat_orientation_sublists for side in sublist]
        #print("Flat genotype after removing all '00' sublists: ", flat_genotype)
    else:
        flat_genotype_sublists = [flat_genotype[i:i + 6] for i in range(0, len(flat_genotype), 6)]
        flat_genotype_sublists = [sublist for sublist in flat_genotype_sublists if not all(side == '00' for side in sublist)]
        flat_genotype = [side for sublist in flat_genotype_sublists for side in sublist]
        flat_orientation = None

    # Simplification #4 --------------------------------------------------------------- 
    # Renumber interfaces respecting their rules of attachment
    new_numbering = {}  # Map old numbers to new ones
    next_number = 1  # Start numbering from 1
    renumbered_flat_genotype = []

    for num in flat_genotype:
        if num == '00':
            renumbered_flat_genotype.append(0)  # Keep '00' as 0
        else:
            if num not in new_numbering:  # If the number is not yet assigned
                new_numbering[num] = next_number  # Assign the next available number
                complement = rules_dict.get(num, None)  # Ensure its complementary partner gets the corresponding pair
                if complement is not None and complement not in new_numbering:
                    new_numbering[complement] = next_number + 1
                next_number += 2  # Increment by 2 to maintain pairing
            renumbered_flat_genotype.append(new_numbering[num])  # Append the renumbered value

    #print("Renumbered genotype: ", renumbered_flat_genotype)
    #print("Renumbered orientation: ", flat_orientation)

    return renumbered_flat_genotype, flat_orientation

#====================================================================================================#
def convert_genostr_to_tile_dict(genotype_string):
    """ 
    Convert string of genotype (list of lists written in a file as string) to tile_dict.

    Args:
        - genotype_string (str): Genotype string read from the all_genotypes file.
    
    Returns:
        - tile_dict (dict): Dictionary with key as the tile name and value as the list of sides.
    """
    import string

    tile_dict = {}
    letters = string.ascii_uppercase  # 'A' to 'Z'

    line = genotype_string
    tiles = line.strip().split('] [')  # Splitting based on the delimiter
    tiles = [tile.replace("[", "").replace("]", "") for tile in tiles]  # Cleaning brackets
    tiles = [tile.split(", ") for tile in tiles]  # Splitting elements inside lists
    tiles = [[elem.strip("'") for elem in tile] for tile in tiles]  # Removing single quotes
    tile_dict = {letters[i]: tiles[i] for i in range(min(len(tiles), len(letters)))}

    return tile_dict

#====================================================================================================#

def convert_genostr_to_genolist(genostr):
    """ 
    Convert genotype string (read from the file) to list of lists.
    [Date: 18 Feb 2025]
    [Author: Prarthana Agrawal]
    [Version: 0.1]

    Args: 
        - genostr (str): Genotype string read from the all_genotypes file.

    Returns:
        - genolist (list): List of lists of genotypes.
    """

    genolist =  [lst.strip("[]").replace("'", "").split(", ") for lst in genostr.split("] [")]

    return genolist


#====================================================================================================#

def flatten_dict(tile_dict):
    """ 
    Flatten the tile_dict to a flat list of sides.

    [Date: 03 Mar 2025]

    Args:
        - tile_dict (dict): Dictionary with key as the tile name and value as the list of sides.

    Returns:
        - flat_list (list): List of all sides of all tiles.
    """

    flat_list = [side for tile in tile_dict.values() for side in tile]

    return flat_list

#====================================================================================================#

def unflatten_list(flat_list):
    """
    Convert a flat list of sides back into a dictionary format with tile names A-Z.

    Args:
        - flat_list (list): The flattened list of sides.

    Returns:
        - tile_dict (dict): Dictionary where keys are tile names (A-Z) and values are lists of sides.
    """
    sides_per_tile = 6
    num_tiles = len(flat_list) // sides_per_tile

    if len(flat_list) % sides_per_tile != 0:
        raise ValueError("Flat list size is not a multiple of 6.")

    tile_names = list(string.ascii_uppercase[:num_tiles])  # Generate 'A', 'B', 'C', etc.
    tile_dict = {tile_names[i]: flat_list[i * sides_per_tile : (i + 1) * sides_per_tile] for i in range(num_tiles)}
    
    # Ensure all values in tile_dict are in two-digit format
    for key in tile_dict:
        tile_dict[key] = [f"{int(side):02d}" for side in tile_dict[key]]
    return tile_dict
 
#====================================================================================================#

def simplify_tiledict(params, tile_dict):
    """ 
    Simplifies the given tile_dict under three simplification rules:
    1. Zero all neutral sides
    2. Zero all sides with no partners
    3. Renumber interfaces respecting their rules of attachment
    #! Note: tiles with all sides '00' are NOT removed here as it may lead to key error.

    Args:
        - params (dict): {'tile_types':tile_types, 'side_types':side_types, 'neutral_sides':neutral_sides, 'self_int_sides':self_int_sides, 'n_rules':n_rules, 'dim':dim}
        - tile_dict (dict): Dictionary with key as the tile name and value as the list of sides.

    Returns:
        - renumbered_tile_dict (dict): Simplified tile_dict.

    """

    from core import rules_func

    # Extract parameters
    neutral_sides = params['neutral_sides']
    
    # Load rules of attachment
    rules_dict = rules_func(params)

    # Flatten the tile_dict to a flat list of sides
    flat_genotype = flatten_dict(tile_dict)
    #print("Given genotype: ", flat_genotype)

    # Simplification #1 ---------------------------------------------------------------
    # Zero all neutral sides
    flat_genotype = ['00' if side in neutral_sides else side for side in flat_genotype]
    #print("All sides after zeroing all neutral sides: ", flat_genotype)

    # Simplification #2 ---------------------------------------------------------------
    # Zero all sides with no partners
    flat_genotype = [side if (side == '00' or rules_dict[side] in flat_genotype) else '00'
    for side in flat_genotype]
    #print("All sides after zeroing all sides with no partners: ", flat_genotype)

    # Simplification #3 --------------------------------------------------------------- 
    # Renumber interfaces respecting their rules of attachment
    new_numbering = {}  # Map old numbers to new ones
    next_number = 1  # Start numbering from 1
    renumbered_flat_genotype = []

    for num in flat_genotype:
        if num == '00':
            renumbered_flat_genotype.append(0)  # Keep '00' as 0
        else:
            if num not in new_numbering:  # If the number is not yet assigned
                new_numbering[num] = next_number  # Assign the next available number
                complement = rules_dict.get(num, None)  # Ensure its complementary partner gets the corresponding pair
                if complement is not None and complement not in new_numbering:
                    new_numbering[complement] = next_number + 1
                next_number += 2  # Increment by 2 to maintain pairing
            renumbered_flat_genotype.append(new_numbering[num])  # Append the renumbered value

    #print("Renumbered genotype: ", renumbered_flat_genotype)
    #print("Renumbered orientation: ", flat_orientation)

    # Convert renumbered flat genotype back to tile_dict
    renumbered_tile_dict = unflatten_list(renumbered_flat_genotype)

    return renumbered_tile_dict

#====================================================================================================#
def get_genotype_for_this_shape(path, filepath, n_tiles, n_sides, tot_splits, target_shape_index, stop_at_first=True):
    """ 
    Get representative genotype for given phenotype/shape.
    [Author: Prarthana Agrawal]
    [Date: 11 Aug 2025]
    [Version: 0.1]

    # Older versions:
        build using get_sample_genotype_indices function but specifically for target shape.

    Args:
        - path (str): path to the folder where data files live.
        - filepath (str): folder/file details within the main directory.
        - n_tiles (int): Number of tile species.
        - n_sides (int): Number of sides/ colors.
        - tot_splits (int): Number of splits for the parallel run.
        - target_shape_index (int): Index of the target shape in the combined valid shapes file.
        - stop_at_first (bool): If True, stop at the first occurrence of the target shape. If False, continue searching for all occurrences.
    
    Returns:
        - list of tuples: Each tuple contains (nsplit, index) where the target shape is found.
    """

    # ==========================================================================================================# 
    # ------------------------- Extract shape from combined valid shapes files -------------------------------- #
    # ==========================================================================================================# 
    """
    Step1: Load combined_valid_shapes data file and extract the final valid shapes.
    """
    with open(path + filepath + 'valid_shapes.txt', 'r') as shape_file:
        for line_index, line in enumerate(shape_file):
            if line_index == target_shape_index:
                coords = eval(line.strip())
                shift_coords = shift_coordinates(coords)
                target_shape = shift_coords  # Get the target shape
                break
        else:
            raise ValueError(f"Line {target_shape_index} not found in file.")

    print("Target shape: ", target_shape)

    # ==========================================================================================================# 
    # -------------------------- Search for target shape in each nsplit_valid_shapes -------------------------- #
    # ==========================================================================================================# 
    """
    Step2: For target shape, find if it exists in each nsplit_valid_shapes. If yes, return its index.
    For example [(0,0,0)] can be shape number #1 in combined_valid_shapes but #3 in nsplit=2.
    This is because in each nsplit, the shapes are numbered in order of appearance.
    """
    nsplit_shape_num_list = list()  # list of tuples (nsplit, shape_num) where the target shape is found
    for nsplit in range(1, tot_splits+1):

        nsplit_filepath = filepath.replace('_combined_', '_nsplit_') + f"{nsplit}_"

        with open(path + nsplit_filepath + 'valid_shapes.txt', 'r') as nsplit_shape_file:
            for shapenum, line in enumerate(nsplit_shape_file, start=1):
                shape2 = eval(line.strip())
                if compare_polycubes(target_shape, shape2) == 1:
                    nsplit_shape_num_list.append((nsplit, shapenum))
                    break

    #===========================================================================================================#
    # Go at this (nsplit, shapenum) in nsplit_shape_type file and get the index at which this shape num exists in shape_type file
    all_nsplit_indices = list()  # list of tuples (nsplit, index) where the target shape is found
    for (nsplit, shapenum) in nsplit_shape_num_list:
        nsplit_filepath = filepath.replace('_combined_', '_nsplit_') + f"{nsplit}_"
        with open(path + nsplit_filepath + 'shape_type.txt', 'r') as nsplit_shape_type_file:
            for index, line in enumerate(nsplit_shape_type_file):
                shape_type = int(line.strip())
                if shape_type == shapenum:
                    all_nsplit_indices.append((nsplit, index))


    return nsplit_shape_num_list, all_nsplit_indices

#===========================================================================================================#

def get_rep_genotype_indices_for_target_complexity(path, filepath, n_tiles, n_sides, tot_splits, complexity_name_for_filter, num_of_groups):
    """ 
    Get representative genotype for each unique phenotype belonging in groups G1 to G5 (simplest 20% to most complex 20% phenotypes).
    This function is similar to get_sample_genotype_indices but specifically for target complexity phenotypes.
    
    [Author: Prarthana Agrawal]
    [Date: 26 Aug 2025]
    [Version: 0.3]

    # Older versions:
        build using get_sample_genotype_indices function but specifically for target shape.

        [Version: 0.1 (11 Aug 2025), 0.2 (25 Aug 2025)]
            - added 'rand' option to randomly pick phenotypes.
            - instead of high/mid/low, use groups G1 to G5 (simplest 20% to most complex 20% phenotypes)

    Args:
        - path (str): path to the folder where data files live.
        - filepath (str): folder/file details within the main directory.
        - n_tiles (int): Number of tile species.
        - n_sides (int): Number of sides/ colors.
        - tot_splits (int): Number of splits for the parallel run.
        - complexity_name_for_filter (str): Name of the complexity file to be used for grouping phenotypes.
        - num_of_groups (int): How many groups to divide the phenotypes into based on complexity. Generally 5.
            1. G1: simplest 20% phenotypes
            2. G2: 20-40% phenotypes
            3. G3: 40-60% phenotypes
            4. G4: 60-80% phenotypes
            5. G5: most complex 20% phenotypes

    Returns:
        - sample_indices_for_all_groups (list): List of sample indices for all groups.
          List of dictionaries sample_indices_dict for each group.
            sample_indices_dict (dict): Dictionary with key as the original shape number in combined file and value as the list of sample indices (tuple of nsplit and index for each phenotype.
            eg {1: [(2,1032), (3,43)], 2: [(8,31)]}
    """

    # ==========================================================================================================# 
    # ---------------------------------- Load combined valid shapes data -------------------------------------- #
    # ==========================================================================================================# 
    """
    Step1: Load combined_valid_shapes data file and extract the final valid shapes.
    """
    shape_file = open(path + filepath + 'valid_shapes.txt', 'r')
    shape_content = shape_file.readlines()

    combined_shapes = [] # final valid shapes from combined file
    for line in shape_content:
        coords = eval(line.strip())
        shift_coords = shift_coordinates(coords) # make coordinates positive
        combined_shapes.append(shift_coords)

    #print("Unique shapes \n", combined_shapes)
    #print("Number of unique shapes: ", len(combined_shapes))

    # Filter shapes based on complexity
    filtered_shapes = []
    filtered_shape_numbers = []

    shape_numbers = np.arange(1, len(combined_shapes) + 1)  # Shape numbers start from 1
    complexity_data = np.loadtxt(path + filepath + complexity_name_for_filter + '.txt', dtype=int)  # Load complexity 
    sorted_complexity, sorted_shape_numbers = sort_together([complexity_data, shape_numbers], reverse=False)  # Sort complexity and shape numbers together in ascending order of complexity values

    num_of_shapes = len(sorted_complexity)

    # Split into groups as evenly as possible
    grouped_complexity = np.array_split(sorted_complexity, num_of_groups)
    grouped_shape_numbers = np.array_split(sorted_shape_numbers, num_of_groups)

    sample_indices_for_all_groups = list()  # list of sample indices for all groups

    for g in range(num_of_groups):
        print(f"Group G{g+1}:")
        print("Complexity range:", min(grouped_complexity[g]), "-", max(grouped_complexity[g]))
        # print("Shape numbers:", grouped_shape_numbers[g])

        filtered_complexity = grouped_complexity[g]
        filtered_shape_numbers = grouped_shape_numbers[g]


        # ==========================================================================================================#
        # ------------For each unique shape, find its corresponding number in each nsplit_valid_shapes ------------ #
        # ==========================================================================================================# 
        """
        Step2: For each filtered unique shape, find its corresponding number in each nsplit_valid_shapes.
        For example [(0,0,0)] can be shape number #1 in combined_valid_shapes but #3 in nsplit=2.
        This is because in each nsplit, the shapes are numbered in order of appearance.
        """

        shape_num_dict = dict() # key: shape number, value: list of corresponding shape numbers in each nsplit
        sample_indices_dict = dict() # key: original shape number in combined file, value: list of sample indices

        for num1 in filtered_shape_numbers:
            shape1 = combined_shapes[num1-1]  # Get the shape corresponding to the shape number

            #print("**************** Shape number: ", num1, "****************")
            match_list = list() # list of corresponding shape numbers in each nsplit

            for nsplit in range(1, tot_splits+1):

                nsplit_filepath = f"{n_tiles}s{n_sides}c_nsplit_{nsplit}_"
                #nsplit_filepath = f"data_files/{n_tiles}s{n_sides}c_nsplit_{nsplit}_"
                nsplit_shape_file = open(path + nsplit_filepath + 'valid_shapes.txt', 'r')
                nsplit_shape_content = nsplit_shape_file.readlines()
                nsplit_shapes = [eval(line.strip()) for line in nsplit_shape_content] # list of shapes in nsplit

                found = 'no'
                for num2, shape2 in enumerate(nsplit_shapes,1):
                    if compare_polycubes(shape1, shape2) == 1: #! make sure the function shifts coordinates first
                        """
                        %matplotlib inline
                        params = dict()
                        fig = plt.figure(figsize=(8, 6))
                        
                        ax1 = fig.add_subplot(121, projection='3d')
                        ax1.view_init(elev=-149, azim=138)
                        plot_all_cubes(params, ax1, shape1, picked_tiles=[], cube_outline='False', axes_lines='False')
                        ax1.set_title(f"Shape {num1}")

                        ax2 = fig.add_subplot(122, projection='3d')
                        ax2.view_init(elev=-149, azim=138)
                        plot_all_cubes(params, ax2, shape2, picked_tiles=[], cube_outline='False', axes_lines='False')
                        ax2.set_title(f"Shape {num2}")
                        
                        plt.show()
                        """
                        found = 'yes'
                        match_list.append(num2) # shape 'type' number in nsplit
                        break

                if found == 'no':
                    match_list.append(0)

            #print(f"Shape {num1} is represented by: \n", match_list)
            shape_num_dict[num1] = match_list
            #print("Match list: ", match_list)

            #----------------------------------------------------------------------------------------------------------#
            # Filter nsplits where the shape is actually present, along with its shape number
            nsplit_matchnum_list = list() # list of tuples (nsplit, shape number in nsplit)
            valid_nsplits = [split for split in range(1, tot_splits+1) if match_list[split-1] != 0]

            for split in valid_nsplits:
                match_num = match_list[split-1]
                nsplit_matchnum_list.append((split, match_num))
            
            #print("(nsplit, shape number in nsplit) \n", nsplit_matchnum_list)

            #==========================================================================================================#    #-----------Search shape_type file and randomly pick indices where desired shape type exists---------------#
            #==========================================================================================================#
            """
            Step3: Randomly pick nsplits from the list of nsplits where the shape is present.
            For each nsplit, randomly pick an index until the shape is found.
            Repeat for sample_size times.

            Step3 (deprecated): Open nsplit_shape_type files and find (nsplit, index) values where the shape type if found.
            we now have a list of nsplits where the shape is present and the corresponding shape number for that nsplit.
            we can then later extract genotypes at this (nsplit, index) values.
            """

            # load shape type files
            #matching_indices = list() # list of tuples (nsplit, index) where the shape number is found
            
            # randomly pick nsplit from the list of nsplits where the shape is present
            random_nsplit_w_shapenum_list = [random.choice(nsplit_matchnum_list)]
            #print("Randomly picked nsplit: ", pick_nsplit_w_shapenum)

            sample_indices = list() # list of tuples (nsplit, index) where the shape number is found
            # open nsplit_shape_type file and randomly pick indices until the shape is found
            for (nsplit, match_num) in random_nsplit_w_shapenum_list:

                #! remove this
                nsplit_filepath = f"{n_tiles}s{n_sides}c_nsplit_{nsplit}_"
                #nsplit_filepath = f"data_files/{n_tiles}s{n_sides}c_nsplit_{nsplit}_"
                #! check this path+
                nsplit_shape_type = open(path + nsplit_filepath + 'shape_type.txt', 'r')

                # Use reservoir sampling to randomly pick indices
                count = 0
                chosen_index = None

                for idx, shape_type in enumerate(nsplit_shape_type,0):
                    if int(shape_type) == match_num:
                        count += 1
                        if random.randint(1, count) == 1:
                            chosen_index = idx

                #! This should not happen if the shape is present in the nsplit_shape_type file
                #! but just in case, raise an error
                if chosen_index is None:
                    raise ValueError(f"Shape number {match_num} not found in nsplit {nsplit} shape_type file.")
            
            sample_indices.append((nsplit, chosen_index))
            # if the required number of samples is more than available, change sample size for this phenotype
            #if sample_size > len(matching_indices):
            #    sample_size = len(matching_indices)
                
            # randomly pick sample_size genotype indices 
            #sample_indices = random.sample(matching_indices, sample_size)
            #print("Sample indices with nsplits: ", sample_indices)

            sample_indices_dict[num1] = sample_indices

        sample_indices_for_all_groups.append(sample_indices_dict)
    return sample_indices_for_all_groups

#===========================================================================================================#
def get_genotype_groups_based_on_complexity(path, filepath, tot_splits, complexity_name_for_filter, num_of_groups=5):
    """ 
    Get groups of genotype (nsplit, index) based on the complexity of their phenotypes.
    #! This is different from get_rep_genotype_indices_for_target_complexity function which gets representative genotype for each unique phenotype in each group. Here we sample genotypes directly instead of sampling phenotypes.

    [Author: Prarthana Agrawal]
    [Date: 17 Sep 2025]
    [Version: 0.1]

    # Older versions:
        build using get_rep_genotype_indices_for_target_complexity function but specifically for sampling genotypes directly.

    Args:
        - path (str): path to the folder where data files live.
        - filepath (str): folder/file details within the main directory.
        - tot_splits (int): Number of splits for the parallel run.
        - complexity_name_for_filter (str): Name of the complexity file to be used for grouping phenotypes.
        - num_of_groups (int): How many groups to divide the genotypes into based on complexity. Generally 5.
            1. G1: simplest 20% phenotypes
            2. G2: 20-40% phenotypes
            3. G3: 40-60% phenotypes
            4. G4: 60-80% phenotypes
            5. G5: most complex 20% phenotypes

    Returns:
        - grouped_rows (list of lists): Outer list = groups, inner list = [(nsplit, index), complexity, complexity_species, lz_complexity]
    """

    # mapping from name -> index in tuple
    col_map = {"complexity": 1, "complexity_species": 2, "lz_complexity": 3}

    if complexity_name_for_filter not in col_map:
        raise ValueError(
            f"Invalid complexity filter: {complexity_name_for_filter}. "
            f"Choose from {list(col_map.keys())}"
        )

    #===========================================================================================================#
    # Collect genotypes (nsplit, index) along with their phenotypic complexities from all nsplit files
    #===========================================================================================================#
    rows = []
    for nsplit in range(1, tot_splits + 1):
        nsplit_filepath = filepath.replace('_combined_', '_nsplit_') + f"{nsplit}_"

        nsplit_shape_type = np.loadtxt(path + nsplit_filepath + 'shape_type.txt', dtype=int)
        nsplit_complexity = np.loadtxt(path + nsplit_filepath + 'complexity.txt', dtype=int)
        nsplit_complexity_species = np.loadtxt(path + nsplit_filepath + 'complexity_species.txt', dtype=int)
        nsplit_lz_complexity = np.loadtxt(path + nsplit_filepath + 'lz_complexity.txt', dtype=int)

        for i, shape_num in enumerate(nsplit_shape_type):
            if shape_num != 0:
                index = shape_num - 1
                comp = nsplit_complexity[index]
                comp_species = nsplit_complexity_species[index]
                lz_comp = nsplit_lz_complexity[index]
                rows.append(((nsplit, i), comp, comp_species, lz_comp))

    #===========================================================================================================#
    # Sort based on the chosen complexity
    #===========================================================================================================#
    sort_index = col_map[complexity_name_for_filter]
    sorted_rows = sorted(rows, key=lambda r: r[sort_index])

    #===========================================================================================================#
    # Split into groups (G1..Gn)
    #===========================================================================================================#
    def split_list(lst, num_groups):
        if num_groups > len(lst):
            num_groups = len(lst)
        n = len(lst)
        group_size = n // num_groups
        remainder = n % num_groups
        groups = []
        start = 0
        for i in range(num_groups):
            end = start + group_size + (1 if i < remainder else 0)
            groups.append(lst[start:end])
            start = end
        return groups

    grouped_rows = split_list(sorted_rows, num_of_groups)

    return grouped_rows


#===========================================================================================================#

def find_genotype_for_this_phenotype(path, n_tiles, n_sides, tot_splits, phenotype_type, target_phenotype, phenotype_name=None, how_many_genotypes=1, save_output=True):
    """ 
    Get genotype(s) for given phenotype from all nsplit files.

    Phenotype maybe: target shape, target size or others.
    #? Currently implemented only shape and size as options.

    Date: 23 Feb 2026
    Author: Prarthana Agrawal
    Version: 0.1

    Older versions:
        v0.0 Built upon genotype_utils/get_sample_genotype_indices

    Args: 
        - path (str): path to the folder where data files live.
        - n_tiles (int): Number of tile species.
        - n_sides (int): Number of sides/ colors.
        - tot_splits (int): Number of splits for the parallel run.
        - dim (int): Dimension of the polycube (2D or 3D).
        - phenotype_type (str): Type of phenotype to search for. Options: 'shape', 'size'.  
        - target_phenotype: The specific phenotype to search for. If phenotype_type='shape', this should be a list of coordinates. If phenotype_type='size', this should be an integer.
        - phenotype_name (str): Name of the phenotype to be used in output file names. Optional, but recommended for clarity. If not provided, a default name based on the phenotype type and target will be used.
        - how_many_genotypes (int/str='all'): Number of genotypes to extract for the target phenotype. Default is 1. 
        - save_output (bool): Whether to save the output in a file. Default is True.
        
        Useful for extracting all genotypes for a specific shape:
            If phenotype type is shape, we could store all genotype indices where the shape is found. 
        Useful for extracting one initial genotype walker for NNSE thermalization:
            If phenotype type is size, we could randomly pick one genotype index from the list of indices, or stop at the first index where shape is found.
    
    Returns:

    """
    #===========================================================================================================#
    #--------------------------------------------- Load packages -----------------------------------------------#
    #===========================================================================================================#
    import numpy as np
    import random
    import linecache
    import string
    from operator import index
    from more_itertools import sort_together

    from utils import shift_coordinates
    from utils.genotype_utils import extract_genotype, extract_orientation
    from symmetry import compare_polycubes
    from plots import plot_all_cubes

    #===========================================================================================================#
    #------------------------------------------- Assertions and checks -----------------------------------------#
    #===========================================================================================================#
    if phenotype_type == 'shape':
        if not isinstance(target_phenotype, list) or not all(isinstance(coord, tuple) for coord in target_phenotype):
            raise ValueError("target_phenotype must be a list of tuples when phenotype_type is 'shape'")
        if phenotype_name is None:
            size = len(target_phenotype)
            phenotype_name = f"shape_{size}"

    elif phenotype_type == 'size':
        if not isinstance(target_phenotype, int):
            raise ValueError("target_phenotype must be an integer when phenotype_type is 'size'")
        if phenotype_name is None:
            phenotype_name = f"size_{target_phenotype}"

    
    #===========================================================================================================#
    #----------------------------------------- Find shape number in each nsplit --------------------------------#
    #===========================================================================================================#
    # First find out what the shape number of the target shape is in the combined valid shapes file. This will be our reference point to find the shape number in each nsplit and then the corresponding genotype indices.
    nsplit_and_indices_dict = {} # to store nsplit and indices where shape is found
    nsplit_numbers_shuffled = list(range(1, tot_splits+1))
    random.shuffle(nsplit_numbers_shuffled) # shuffle the nsplit numbers to ensure random sampling

    for nsplit in nsplit_numbers_shuffled:
        
        filepath = f"data_files/{n_tiles}s{n_sides}c_nsplit_{nsplit}_"
        nsplit_shape_file = open(path + filepath + 'valid_shapes.txt', 'r')
        nsplit_shape_content = nsplit_shape_file.readlines()
        nsplit_shapes = [eval(line.strip()) for line in nsplit_shape_content] # list of shapes in nsplit

        found = 'no'
        for num2, shape2 in enumerate(nsplit_shapes,1):

            if phenotype_type == 'shape':
                target_shape = shift_coordinates(target_phenotype) # shift to origin
                shape2 = shift_coordinates(shape2) # shift to origin

                if compare_polycubes(target_shape, shape2) == 1: 
                    found = 'yes'
                    target_shape_number = num2
                    break

            # here if the phenotype type is size, we look for first matching size in the list of all shapes. This shape number is enough to get the genotype index from shape type file. We can later modify this to look for all shapes with the target size and then get all their shape numbers and corresponding genotype indices.
            elif phenotype_type == 'size':
                
                target_size = target_phenotype
                shape2_size = len(shape2)

                if target_size == shape2_size:
                    found = 'yes'
                    target_shape_number = num2
                    break

            else: 
                raise ValueError("Invalid phenotype type. Choose from 'shape' or 'size'.")

        if found == 'no':
            continue

        #===========================================================================================================#
        #------------------- Search shape_type file and pick indices where desired phenotye exists -----------------#
        #===========================================================================================================#
        # Now, we know the shape number of the target phenotype (shape or size) in the nsplit. We can now go to the shape_type file and find indices where this shape number is found. We can pick one of these indices for nnse thermalization purpose. 
        nsplit_genotype_indices = []
        nsplit_shape_type_file = np.loadtxt(path + filepath + 'shape_type.txt')
        for idx, shape_num in enumerate(nsplit_shape_type_file):
            if int(shape_num) == target_shape_number:
                nsplit_genotype_indices.append(idx)
                # If you only want one genotype per phenotype, stop here. If you want more, keep going until you find the required number of genotypes or exhaust the shape_type file.
                if how_many_genotypes != 'all' and len(nsplit_genotype_indices) >= how_many_genotypes:
                    break

        nsplit_and_indices_dict[nsplit] = nsplit_genotype_indices

        #! this needs to be checked further
        #! here we go through all nsplits, and then pick required number of genotypes from that nsplit. but if we want it to be randomly sampled (nsplit, index) across all nsplits, this needs to be modified. we can first collect all (nsplit, index) where the shape is found across all nsplits, and then randomly pick required number of genotypes from this combined list.
        if how_many_genotypes != 'all' and len(nsplit_and_indices_dict) >= how_many_genotypes:
            break

    #print("nsplit and genotype indices where shape is found:")
    #print(nsplit_and_indices_dict)

    #print("total number of genotypes found for the shape across all nsplits:")
    #total_genotypes = sum(len(indices) for indices in nsplit_and_indices_dict.values())
    #print(total_genotypes)  

    #===========================================================================================================#
    #------------------------------------------ Saving to output files -----------------------------------------#
    #===========================================================================================================#
    # Save the nsplit and its genotype indices in a file
    if save_output:
        output_file = f"{phenotype_name}_genotype_indices.txt"
        with open(output_file, 'w') as f:
            for nsplit, indices in nsplit_and_indices_dict.items():
                f.write(f"nsplit: {nsplit}, indices: {indices}\n")

        print(f"nsplit and genotype indices saved to {output_file}")

    # Extract genotypes at these indices
    output_file = f"{phenotype_name}_extracted_genotypes.txt"
    output_file2 = f"{phenotype_name}_extracted_orientations.txt"

    # Clear the output files if they exist
    if save_output:
        open(output_file, "w").close()
        open(output_file2, "w").close()

    genotype_list = []
    orientation_list = []

    for nsplit, genotype_indices in nsplit_and_indices_dict.items():
        filepath = path+f"data_files/{n_tiles}s{n_sides}c_nsplit_{nsplit}_"

        for index in genotype_indices:
        
            genotype = extract_genotype(filepath, index)
            orientation = extract_orientation(filepath, index)

            if save_output:
                with open(output_file, 'a') as f1:
                    f1.write(genotype + '\n')

                with open(output_file2, 'a') as f2:
                    f2.write(str(orientation) + '\n')

            else:
                genotype_list.append(genotype)
                orientation_list.append(orientation)

    if not save_output:
        return genotype_list, orientation_list
    else:
        return f1.name, f2.name