"""
# Description: Contains functions to check for genotype validity.
    - A genotype is considered 'valid' only when it gives rise to 'bounded' and 'deterministic' outputs.

# Functions:
    - valid_sol_checker(params, assembly_settings, tile_dict, orient_dict, sol_stats): to check if the rule gives a bounded, deterministic shape.

# Dependencies:
    - copy
    - matplotlib.pyplot
    - mpl_toolkits.mplot3d.Axes3D
    - utils.return_lengths
    - core.assembly_func
    - symmetry.compare_polycubes 
    - plots.plot_all_cubes
"""

# Import packages
import copy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from utils import return_lengths
from core.self_assembly import assembly_func
from symmetry import compare_polycubes
from plots import plot_all_cubes


def valid_sol_checker(params, assembly_settings, tile_dict, orient_dict, sol_stats, return_all_complexities=False):
    """  
    Checks if the rule gives a bounded, deterministic shape.
    Returns 1 if solution is valid else returns 0, along with the tile coordinates

    Version: 0.4 (created on 01 Oct 2025)

    Old Versions:
        # 0.3 (24 July 2025)
            - added return_all_complexities flag to return complexity from all k-runs instead of just minimum complexity. default is False to maintain backward compatibility.
             if return_all_complexities is True, then complexity_list is returned instead of min(complexity_list). This is to analyze the variation in complexity across different k-runs for the same solution

        # 0.2 (created on 04 Feb 2025)
            - correction: corrected return of tile coordinates and picked tiles
            there was a mismatch in tile_coord from last k-run and picked_tiles from first k-run
            now, both are from first k-run

        # 0.1 (created in Jan 2025)
            - modification: enabled unseeded assembly

    Args:
        - params (dict): {'tile_types':tile_types, 'side_types':side_types, 'neutral_sides':neutral_sides, 'self_int_sides':self_int_sides, 'n_rules':n_rules, 'dim':dim}
            tile_types : number of tile species
            side_types : number of side types
            neutral_sides : number of neutral sides
            self_int_sides : number of self-interacting sides
            n_rules : number of rules
            dim : dimension of the polycube assembly
        - assembly_settings (dict): {'assembly_type': assembly_type, 'Dmax': Dmax, 'max_tiles': max_tiles, 'kmax': kmax}
            assembly_type : type of assembly ('unseeded' vs 'seeded)
            Dmax : linear threshold
            max_tiles : number of tiles threshold
            kmax : loop size for non-determinism check
        - tile_dict (dict): dictionary of all tile species with their interfaces.
        - orient_dict (dict): dictionary of all tile orientations.
        - sol_stats (dict): {'UBD': UBD, 'ND': ND, 'valid': valid}
            UBD : number of unbounded solutions
            ND : number of non-deterministic solutions
            valid : number of valid solutions
        - return_all_complexities (bool): if True, returns all complexities of the valid solution instead of just the minimum complexity.

    Returns:
        - output (int): 1 if solution is valid else 0
        - tile_coord (list): coordinates of all placed tiles
        - picked_tiles (list): list of all picked tiles
        - min(complexity_list) (int): minimum complexity of the valid solution
        or complexity_list (list): list of all complexities of the valid solution if return_all_complexities is True
        - sol_stats (dict): {'UBD': UBD, 'ND': ND, 'valid': valid}

    Raises:
        - Exception: If new choice line is not possible.

    """
    
    Dmax = assembly_settings['Dmax']
    max_tiles = assembly_settings['max_tiles']
    kmax = assembly_settings['kmax']
    assembly_type = assembly_settings['assembly_type']

    first_choice_line = []

    #? ========================================== Unseeded Assembly ================================================ #
    if assembly_type == 'unseeded':
        tile_coord, complexity, _, _, picked_tiles = assembly_func(params, assembly_settings, tile_dict, orient_dict, first_choice_line) # the _ are lineage and choice_tree which are not needed here
        complexity_list = [complexity] # stores the complexity from each k-run
        
        (len_x, len_y, len_z) = return_lengths(tile_coord)[0] # width and height of assembly
        
        # ------------------------------------ check for unboundedness ---------------------------------- #
        if len_x >= Dmax or len_y >= Dmax or len_z >= Dmax or len(tile_coord) >= max_tiles: 
            sol_stats['UBD'] = sol_stats['UBD'] + 1; output = 0
            #return output, None, None, None #! newly uncommented
            return output, tile_coord, picked_tiles, min(complexity_list), sol_stats # returns tile coord and picked tiles from the very first k-run
            #! does min(complexity_list) serve any purpose here?

        # ------------------------------------ check for determinism ---------------------------------- #
        else:
    
            very_first_picked_tiles = picked_tiles # to return as output of the function (first k-run picked tiles)
            very_first_shape = tile_coord # to return as output of the function (first k-run shape)
            first_shape = tile_coord # to compare with subsequent k-runs

            count = 1
            for _ in range(kmax):
                choice_line = []
                tile_coord, complexity, _, _, picked_tiles = assembly_func(params, assembly_settings, tile_dict, orient_dict, choice_line)

                complexity_list.append(complexity)
                second_shape = tile_coord.copy()

                if compare_polycubes(second_shape, first_shape) == 1:  # if two shapes are equivalent under symmetry operations  
                    count = count + 1
                    first_shape = second_shape.copy()
                    second_shape = [] # to prevent garbage values filling in for the next iteration
                else:
                    sol_stats['ND'] = sol_stats['ND'] + 1; output = 0
                    #return output, None, None, None # fisrt_shape, min(complexity_list)
                    if return_all_complexities:
                        return output, very_first_shape, very_first_picked_tiles, complexity_list, sol_stats
                    else:
                        return output, very_first_shape, very_first_picked_tiles, min(complexity_list), sol_stats
                        #break # non-deterministic

            if count == kmax + 1:
                output = 1
                sol_stats['valid'] = sol_stats['valid'] + 1
                if return_all_complexities:
                   return output, very_first_shape, very_first_picked_tiles, complexity_list, sol_stats 
                else:
                    return output, very_first_shape, very_first_picked_tiles, min(complexity_list), sol_stats
                #! Correction 1 oct 2025: previously returned first_picked_tiles
                #! this led to prints tile coord from last k-run while picked tiles from first run
            else:
                raise Exception("Error: ND shape should have exited loop earlier. Something's wrong!")

    #? ========================================= Seeded Assembly ================================================ #
    tile_coord, complexity, lineage, first_tree, picked_tiles = assembly_func(params, assembly_settings, tile_dict, orient_dict, first_choice_line)
    full_tree = copy.deepcopy(first_tree) # shallow copy will affect value of full_tree if first_tree is sliced
    terminals = [list(lineage)] # lineage returned is a fully explored path aka terminal end
    #p#print("Terminals {} \nFirst Tree {} \nLineage {}".format(terminals, first_tree, lineage))
      
    complexity_list = [complexity]

    # To plot the polycube 
    #fig = plt.figure()
    #ax1 = fig.add_subplot(111, projection='3d')
    #ax1.view_init(elev=-149, azim=138)
    #plot_all_cubes(params, ax1, tile_coord, picked_tiles, cube_outline='False', axes_lines='False')

            
    # ============================ Check for boundedness ========================== #    
    #* position of this has been modified. earlier len(first_tree) < 2 was tested before this-- misclassified UBD sols
    (len_x, len_y, len_z) = return_lengths(tile_coord)[0] # width and height of assembly
    if len_x >= Dmax or len_y >= Dmax or len_z >= Dmax or len(tile_coord) >= max_tiles: 
        sol_stats['UBD'] = sol_stats['UBD'] + 1; output = 0
        #return output, None, None, None #! newly uncommented
        if return_all_complexities:
            return output, tile_coord, picked_tiles, complexity_list, sol_stats
        else:
            #! does min(complexity_list) serve any purpose here?
            return output, tile_coord, picked_tiles, min(complexity_list), sol_stats 

    else: # ======================== Check for determinism ========================= #
        if len(first_tree) < 2: # len = 0 means seed tile had no attachments, 1 means no alternative choices are there
            sol_stats['valid'] = sol_stats['valid'] + 1; output = 1
            #return output, None, None, None #! newly uncommented
            if return_all_complexities:
                return output, tile_coord, picked_tiles, complexity_list, sol_stats
            else:
                return output, tile_coord, picked_tiles, complexity, sol_stats  #! PLEASE INSPECT IF THIS LINE IS CORRECT #! complexity vs min(complexity_list)

        filtered_nodes = list(filter(lambda node: node not in terminals, first_tree))
        choice_line = filtered_nodes[0] #random.choice(filtered_nodes) 
        #p#print("First Choice Line=", choice_line)

        first_shape = tile_coord
        very_first_picked_tiles = picked_tiles # to return as output of the function (first k-run picked tiles)
        very_first_shape = tile_coord # to return as output of the function (first k-run shape)
        count = 1

        while True: # loop over min(all possible, kmax) paths
            tile_coord, complexity, lineage, choice_tree, picked_tiles = assembly_func(params, assembly_settings,tile_dict, orient_dict, choice_line)
            choice_line = []
            terminals.append(lineage)
            #p#print("Terminals {} \nChoice Tree {} \nLineage {}".format(terminals, choice_tree, lineage))
            new_values = list(filter(lambda val: (val not in full_tree), choice_tree)) # New paths revealed
            #p#print("New Values not in full tree =", new_values)

            #? Filter new paths which are not present in full_tree already (even as subsets)
            remove = list()
            for val in new_values:
                for old in full_tree:
                    if val == old[:len(val)]: #if new value is subset of existing value in fulltree, remove it
                        remove.append(val)
                        break
            
            for val in remove:
                new_values.remove(val)
            #p#print("New Values =", new_values)
            #p#print("Full tree", full_tree)

            full_tree.extend(new_values) # Add "truly" new paths to full_tree
            #p#print("Full tree after adding new paths", full_tree)

            #? Filter out values in full_tree which are subset of one another (these paths are subset of another path)
            data = full_tree
            filtered_data = []
            for i in range(len(data)):
                item = ','.join(data[i]) # create a combined string
                present_in_others = False
                for j in range(len(data)):
                    other_item = ','.join(data[j])  # create a combined string
                    if i != j: # distinct elements
                        if (item[-1] == '_' and item == other_item[:len(item)]) or  (item[-1] != '_' and item+',' == (other_item+',')[:len(item)+1]): # ['5_1'] is not subset of ['5_10']; ['1_'] is subset of ['1_4']
                            present_in_others = True
                            break
                if not present_in_others:
                    filtered_data.append(data[i])
            full_tree = filtered_data
            #p#print("Full tree after removing redundant values", full_tree)

            #============================================== SHAPE COMPARISON ==========================================#
            complexity_list.append(complexity)
            second_shape = tile_coord.copy()

            #fig = plt.figure()
            #ax1 = fig.add_subplot(111, projection='3d')
            #ax1.view_init(elev=-149, azim=138)
            #plot_all_cubes(params, ax1, tile_coord, picked_tiles, cube_outline='False', axes_lines='False')

            if compare_polycubes(second_shape, first_shape) == 1:  # if two shapes are equivalent under symmetry operations  
                count = count + 1
                first_shape = second_shape.copy()
                second_shape = [] # to prevent garbage values filling in for the next iteration
            else:
                sol_stats['ND'] = sol_stats['ND'] + 1; output = 0
                if return_all_complexities:
                    return output, very_first_shape, very_first_picked_tiles, complexity_list, sol_stats
                else:
                    #return output, None, None, None # fisrt_shape, min(complexity_list)
                    return output, very_first_shape, very_first_picked_tiles, min(complexity_list), sol_stats
                #break # non-deterministic

            #? if all elements in full tree are fully explored  or reached kmax iterations then end
            if (len(full_tree) < kmax and sorted(terminals) == sorted(full_tree)) or  (len(full_tree) >= kmax and len(terminals) == kmax):
                #p#print_msg_box("ALL PATHS EXPLORED")
                #p#print("Full Tree", len(full_tree), "\n", full_tree)
                output = 1; sol_stats['valid'] = sol_stats['valid'] + 1
                #p#print("len(terminals)", len(terminals), "len(full_tree)", len(full_tree))
                if return_all_complexities:
                    return output, very_first_shape, very_first_picked_tiles, complexity_list, sol_stats
                else:
                    return output, very_first_shape, very_first_picked_tiles, min(complexity_list), sol_stats
            else:         
                #? Choose a new choice line for next iteration
                filtered_nodes_deepcopy = [copy.deepcopy(node) for node in full_tree if node not in terminals]
                if len(filtered_nodes_deepcopy) <= 0:
                    print("Tile_dict", tile_dict)
                    raise Exception("new choice line not possible")
                
                choice_line = filtered_nodes_deepcopy[0] # random.choice(filtered_nodes_deepcopy)
                #p#print("Next Choice Line=", choice_line)
