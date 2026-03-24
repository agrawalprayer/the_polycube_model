"""
# Description: Contains the main code for simulating self-assembly of polycubes.

    - Number of tile species (Nt) and number of colors (c) are user-defined parameters.
    - The assembly starts with a seed tile and proceeds by attaching tiles to the polycube.
    - The assembly process is terminated if the size of the polycube exceeds a threshold (Dmax) or the number of tiles exceeds a threshold (max_tiles).
    - Only valid solutions are considered, i.e., those which are bounded and deterministic.

# Functions:
    - assembly_func(params, assembly_settings, tile_dict, orient_dict, choice_line): simulates self-assembly of tiles.

# Dependencies:
    - random
    - more_itertools
    - core.rules_func 
    - core.zero_sides
    - core.origin_finder
    - core.path_finder
    - utils.extract_underscore
    - utils.return_lengths   
"""

# Import packages
import random
from more_itertools import sort_together
from core import rules_func, zero_sides, origin_finder, path_finder
from utils import extract_underscore, return_lengths

def assembly_func(params, assembly_settings, tile_dict, orient_dict, choice_line):
    """ 
    Self-assembly of polyominoes.

    Version: 0.3 (created on 03 Mar 2025)

    Old Versions:

        # 0.2 (created on 04 Feb 2025)
            - modification: Zero all sides with no complementary partners

        # 0.1 (created in Jan 2025)
            - modification: enabled unseeded assembly

    Args:
        - params (dict): {'tile_types':tile_types, 'side_types':side_types, 'neutral_sides':neutral_sides, 'self_int_sides':self_int_sides, 'n_rules':n_rules, 'dim':dim}
            tile_types (list): List of tile names eg ['A', 'B'].
            side_types (list): List of interface types eg ['00', '01', '02'].
            neutral_sides (list): List of neutral sides eg ['00', '03'].
            self_int_sides (list): List of self-interacting sides eg ['03'].
            n_rules (int): Number of rules to generate.
            dim (int): dimensionality of the assembly process (2 or 3)
        - assembly_settings (dict): {'assembly_type':assembly_type, 'Dmax':Dmax, 'max_tiles':max_tiles, 'kmax':kmax}
            assembly_type: 'seeded' vs 'unseeded'
            Dmax: linear threshold size of polycube growth
            max_tiles: threshold of maximum number of tiles for growth termination
            kmax: number of checks for non-determinism
        - tile_dict (dict): dictionary mapping tile names with interfaces.
        - orient_dict (dict): dictionary mapping tile names with orientations.
        - choice_line (str): eg indicate the choices of t1 and t2 to follow in the assembly eg'23_4'.

    Returns:
        - tile_coord (list of tuples): list of tile-origins of polyomino
        - complexity (int): number of used interfaces
        - lineage (list): Stores the current choices line/path. Choice of open sides of tile1 and complementary tile2.
        - choice_tree (list): A tree of all turning points and choices made in this iteration.
        - picked_tiles (list): list of names of tiles placed in the polyomino

    Raises:
        - ValueError: Invalid value of dimension.
        """
    # Unpack parameters
    tile_types = params['tile_types']
    neutral_sides = params['neutral_sides']
    dim = params['dim']

    assembly_type = assembly_settings['assembly_type']
    Dmax = assembly_settings['Dmax']
    max_tiles = assembly_settings['max_tiles']

    if dim not in [2,3]: raise ValueError("Invalid value of dim")

    rules_dict = rules_func(params) # rules of attachment
    index_rules_dict = dict({0:2,1:3,2:0,3:1,4:5,5:4}) # tile1's index 0 can attach with tile2's index 2 and so on..

    #?================================================== Seed tile ================================================#
    if assembly_type == 'seeded': seed_tile = tile_types[0] # Using seeded tile assembly i.e. first tile as seed tile
    elif assembly_type == 'unseeded': seed_tile = random.sample(tile_types,1)[0] # [0] here extracts list element
    #! CANNOT USE PATH_FINDER FOR UNSEEDED ASSEMBLY as no seedtile choices are recorded 
    #p#print("Seed tile is", seed_tile, tile_dict[seed_tile])
    picked_tiles = [seed_tile]

    #?========================= Store all available "placed" sides open to further attachment ======================#
    available_sides = list(tile_dict[seed_tile].copy()) # list of open sides # copy is created to avoid overwriting
    available_orient = orient_dict[seed_tile] # string of orientations of open sides
    origin0 = (0,0,0); tile_coord = [origin0] # list of tile origins

    used_interfaces = set(available_sides.copy()) # Used interfaces list to find complexity 
    #used_interfaces = sorted(set(available_sides.copy()), key=available_sides.copy().index) #! sorting set to preserve index--check if needed!!
    #p#print("Open sides = {}, \nTile Coordinates = {}".format(available_sides,tile_coord))
    #p#print("Orientation of open sides", available_orient)

    #?================================= Store all sides (pool of tile2-sides) =====================================#
    all_sides = list()
    for tile in tile_types:
        all_sides.extend(list(tile_dict[tile]))

    #p#print("All tile-2 sides before zeroing =", all_sides)
    all_sides = ['00' if side in neutral_sides else side for side in all_sides]
    #! NEW LINE ADDED - ZEROES SIDES WITH NO COMPLEMENTARY SIDE IN GIVEN GENOTYPE --  REQUIRES TESTING
    all_sides = [side if (side == '00' or rules_dict[side] in all_sides) else '00' for side in all_sides]
    #p#print("All tile-2 sides after zeroing =", all_sides)

    #?========================= Extract indices of non-neutral sides from seed tile ===============================#
    relevant_indices = list(filter(lambda pos:(available_sides[pos] not in neutral_sides and rules_dict[available_sides[pos]] in all_sides),range(len(available_sides)))) # remove open sides which do not have comp side in tile2allsides
    #p#print("Relevant Indices {}".format(relevant_indices))
    
    #?======================== If seed_tile cannot attach to anything, terminate the function =====================#
    if len(relevant_indices) == 0:
        complexity = len(used_interfaces)
        return tile_coord, complexity, [], [], picked_tiles #! What should be the value of lineage and choice tree?

    #?============================================== Assembly =====================================================#
    n_open_pos = 1; #! Should this be initialized as n_open_pos = len(relevant_indices)?
    cutoff = False; first_run = True
    while n_open_pos > 0 and cutoff == False: # run till no open sides left or size exceeds threshold D_max
        #p#print("="*50)
        #p#print("Choice line", choice_line)

        if assembly_type == 'unseeded': choice_line = [] # no choice tree for unseeded assembly #! new addition

        # ================================================ First Run ============================================ #
        if first_run == True:
            #p#print_msg_box("First run")
            # -------------------------------------------- Tile1 ------------------------------------------------ #
            # if a pre-determined 'choice_line' is passed, use that to dictate choice of open sides else pick randoml
            if len(choice_line)==0: random_index = random.sample(relevant_indices,1)[0] #relevant_indices[0]
            else: 
                random_index = extract_underscore(choice_line[0], 'first') 
            random_side = available_sides[random_index] # return side at random index (tile1)
            t1_orient = available_orient[random_index] # return orientation of the chosen random side
            #p#print("Random index = {}, Random side = {}".format(random_index, random_side)) 
            #p#print_msg_box("random index and orientation = "+str(random_index)+" "+str(t1_orient))

            # --------------------------------------------- Tile2 ------------------------------------------------ #
            comp_side = rules_dict[random_side] # return complementary side of the random side
            # indices where complementary side is present in all_sides
            comp_tile2_indices = list(filter(lambda pos:(all_sides[pos]==comp_side), range(0, len(all_sides))))
            
            if len(choice_line)==0 or choice_line[0][-1] == "_": # if no choice line is given for choosing tile2
                t2index = random.sample(comp_tile2_indices,1)[0] #comp_tile2_indices[0]    
            else: 
                t2index = extract_underscore(choice_line[0], 'second') # if choice line is given
            if len(choice_line) > 0: del choice_line[0]  
            #p#print("Complementary side = {} \nComplementary tile2 indices = {} \nChosen t2 index = {}".format(comp_side, comp_tile2_indices, t2index))
            #p#print_msg_box("t2 index = "+str(t2index))
            
            # ------------------------------------------- Choice Tree --------------------------------------------- #
            lineage = [str(random_index)+"_"+str(t2index)] # stores current choices made
            other_choices = list(filter(lambda ind: (ind != random_index), relevant_indices)) # other open side choices
            # dynamically create a tree with new branches for other choices
            if len(other_choices) > 0: choice_tree = [[str(oc)+"_"] for oc in other_choices]
            else: choice_tree = []
            for t2in in comp_tile2_indices:
                choice_tree.append([str(random_index)+"_"+str(t2in)]) # Open side with possible t2 index choices
            #p#print("Choice tree {}, lineage {}".format(choice_tree, lineage))
            first_run = False

        # ================================================== Next Runs ============================================ #
        else:

            if assembly_type == 'unseeded':

                random_index = random.sample(relevant_indices,1)[0] #relevant_indices[0] 
                random_side = available_sides[random_index] # return side at random index (tile1)
                t1_orient = available_orient[random_index] # return orientation of the chosen random side
                #! new line added
                
                comp_side = rules_dict[random_side] # return complementary side of the random side
                comp_tile2_indices = list(filter(lambda pos:(all_sides[pos]==comp_side), range(0, len(all_sides)))) #allcompindices
                t2index = random.sample(comp_tile2_indices,1)[0] #comp_tile2_indices[0] 
                
                lineage = []
                choice_line = []
                choice_tree = []

            else:

                parent = lineage.copy()
                random_index, t1_orient, t2index, lineage, choice_line, choice_tree = path_finder(available_sides, available_orient, all_sides, relevant_indices, parent, lineage, choice_line, choice_tree, t2index, comp_tile2_indices, rules_dict)
        
        # If there are no complementary sides for this chosen open side - ideally this should not happen!!
        if len(comp_tile2_indices) == 0:
            raise Exception("No Complementary side found for this open side")

        # =============================================== ACTUAL ASSEMBLY PROCESS =================================== #
        d = 6 #2*dim deprecated # 2D assemblies are also written as 3D assemblies but with last two indices empty
        tile2_num = int(t2index/d) # return tile number of corresponding chosen index
        tile2 = chr(tile2_num+65) # return tile name of corresponding chosen index
        picked_tiles.append(tile2)
        tile2_interface = tile_dict[tile2] # interface of tile2
        tile2_orients = orient_dict[tile2] # orientations of interfaces of tile2
        tile2_index = t2index % d # relative index of comp side on tile2_interface {0,1,2,3,4,5}
        #p#print("t2 index = {}, tile2_num = {}, tile2 = {}, tile2_interface = {}, tile2_index = {}".format(t2index, tile2_num, tile2, tile2_interface, tile2_index))

        tile2_exp_index = index_rules_dict[random_index % d] # expected index of comp side for attachment
        diff = tile2_exp_index - tile2_index # difference between expected and actual index of comp side
        rotated_interface = tile2_interface.copy() # to store the rotated set of interfaces
        #print("Actual index={}, Expected index={}, Diff={}".format(tile2_index, tile2_exp_index, diff))
        used_interfaces.update(tile2_interface) # store all unique interfaces used --to calculate complexity

        # list of all rotations of a cube with respective orientations
        rot_orient_dict = dict({
            '123456':['OUIDLR','LIROUD','IDOURL','ROLIDU'], '234156':['LIRODU','IUODLR','ROLIUD','ODIURL'], 
            '341256':['IUODRL','LORIDU','ODIULR','RILOUD'], '412356':['LORIUD','OUIDRL','RILODU','IDOULR'],
            '153642':['LIRODU','IUODLR','ROLIUD','ODIURL'], '536142':['IUODRL','LORIDU','ODIULR','RILOUD'], 
            '361542':['LORIUD','OUIDRL','RILODU','IDOULR'], '615342':['OUIDLR','LIROUD','IDOURL','ROLIDU'],
            '143265':['IDOULR','RILODU','OUIDRL','LORIUD'], '432165':['ROLIDU','IDOURL','LIROUD','OUIDLR'], 
            '321465':['ODIURL','ROLIUD','IUODLR','LIRODU'], '214365':['RILOUD','ODIULR','LORIDU','IUODRL'],
            '163524':['RILOUD','ODIULR','LORIDU','IUODRL'], '635124':['IDOULR','RILODU','OUIDRL','LORIUD'], 
            '351624':['ROLIDU','IDOURL','LIROUD','OUIDLR'], '516324':['ODIURL','ROLIUD','IUODLR','LIRODU'],
            '254613':['IDOULR','RILODU','OUIDRL','LORIUD'], '546213':['ROLIDU','IDOURL','LIROUD','OUIDLR'], 
            '462513':['ODIURL','ROLIUD','IUODLR','LIRODU'], '625413':['RILOUD','ODIULR','LORIDU','IUODRL'],
            '264531':['OUIDLR','LIROUD','IDOURL','ROLIDU'], '645231':['LIRODU','IUODLR','ROLIUD','ODIURL'], 
            '452631':['IUODRL','LORIDU','ODIULR','RILOUD'], '526431':['LORIUD','OUIDRL','RILODU','IDOULR']
        })

        # Filter which rotation of interfaces is necessary to attach tile2
        relevant_rotations = [rot for rot in rot_orient_dict.keys() if rot[tile2_exp_index] == str(int(tile2_index)+1)]
        comp_side_orient = tile2_orients[tile2_index] # orientation of tile2 index which is to be bonded
        starting_orients = ['OUIDLR','LIROUD','IDOURL','ROLIDU'] # starting orientation of any cube {0,1,2,3}
        init_cube_orient = [orient for orient in starting_orients if orient[tile2_index] == comp_side_orient][0]
        orient_num = starting_orients.index(init_cube_orient) # find which orientation number for aligning tile2index
        rot = [rot for rot in relevant_rotations if rot_orient_dict[rot][orient_num][tile2_exp_index] == t1_orient][0] 
        # choose that rotation and orientation for which tile2 index orients as required eg at index 3 'Up'
        rotated_interface = [tile2_interface[int(rot[i])-1] for i in range(d)] 
        #p#print("Relevant Rotations {} \nChosen Rotation {} \nRotated Interface {}".format(relevant_rotations, rot, rotated_interface))
    
        update_orient_temp = '' #temporary variable to store new orientations
        for oi in range(d):
            orient = tile2_orients[oi] #orientation of tile2 interface at given index
            or_num = starting_orients.index([ele for ele in starting_orients if ele[oi]==orient][0]) #which orientation rule matches with orientation of present tile2 interface
            rot_sort, ori_sort = sort_together([rot, rot_orient_dict[rot][or_num]], reverse=False) #rewrite the orientation and rotation in ascending order so as to match indices
            new_orient = ori_sort[oi] #under chosen orientation rule, the new orientation should be this
            update_orient_temp = update_orient_temp + new_orient #store new orientations in ascending order
        updated_orient = [update_orient_temp[int(rot[i])-1] for i in range(d)] #reorder the new orients as per the rotation rule
        updated_orient = ''.join(updated_orient) # create a concatenated string from all of the orientations
        #print(updated_orient)

        # Zero the sides which are: attached, neutral, blocked (because of position)
        rotated_interface[tile2_exp_index] = '00' # attached side (tile2)
        available_sides[random_index] = '00' # attached side (tile1)

        # If a side is blocked, rewrite its side as '00'
        origin = tile_coord[int(random_index/d)]  # tile1-origin (6 interfaces correspond to one tile and one origin)
        t2_origin = origin_finder(origin,random_index%d, params) # tile2-origin

        if t2_origin in tile_coord: 
            raise Exception("Error in the code: This position is blocked!") #! CHECK IF NECESSARY
        rotated_interface, available_sides = zero_sides(t2_origin, tile_coord, rotated_interface, available_sides) # zero indirectly blocked neighbouring tiles

        tile_coord.append(t2_origin) # Store origin of attached tile
        available_sides.extend(rotated_interface) # Add the new tile's sides in open sides
        available_orient = available_orient + updated_orient # Add the new orientations of tile2 to available_orient
        #p#print("Available sides = {}, \nTile Coordinates = {}".format(available_sides,tile_coord))
        #p#print("Orientation of open sides", available_orient)

        # Zero the neutral sides
        available_sides = ['00' if x in neutral_sides else x for x in available_sides]
        #p#print("Available sides after zeroing = {}".format(available_sides))

        # Non-zero sides with available partners are relevant sides
        relevant_indices = list(filter(lambda pos:(available_sides[pos] not in neutral_sides and rules_dict[available_sides[pos]] in all_sides),range(len(available_sides)))) # remove open sides which do not have comp side in tile2allsides
        n_open_pos = len(relevant_indices) # Number of open sides
        #p#print("Relevant indices {} \nn_open_pos {}".format(relevant_indices, n_open_pos))
        
        # =================================== Criteria for unbounded assembly ======================================== #
        len_x, len_y, len_z = return_lengths(tile_coord)[0]
        cutoff = any(li >= Dmax for li in [len_x, len_y, len_z]) or len(picked_tiles) > max_tiles

    complexity = len(used_interfaces)
    return tile_coord, complexity, lineage, choice_tree, picked_tiles
