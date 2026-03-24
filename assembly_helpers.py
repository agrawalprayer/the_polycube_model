"""
# Description: Contains functions "helping" in simulating self-assembly process.

# Functions:
    - origin_finder(t1_origin, index, params): Determines tile2's origin relative to tile1's origin.
    - zero_sides(t2_origin, tile_coord, rot_interfaces, available_sides): Zeroes the touching faces of indirectly-blocked neighbouring tiles.
    - path_finder(available_sides, available_orient, all_sides, relevant_indices, parent, lineage, choice_line, choice_tree, t2index, comp_tile2_indices, rules_dict): Generates a choice tree while tracing a particular path (lineage).

# Dependencies:
    - random (path_finder)
    - utils.extract_underscore (path_finder)
"""

#==========================================================================================================#
def origin_finder(t1_origin, index, params):
    """
    Determines tile2's origin relative to tile1's origin based on a positional index (012345 for NESWBF).
    tile1 is the tile already placed. tile2 is the incoming tile for attachment.

    Args:
        - t1_origin (tuple): origin of tile1.
        - index (int): positional-index 012345:NESWBF.
        - params (dict): {'tile_types':tile_types, 'side_types':side_types, 'neutral_sides':neutral_sides, 'self_int_sides':self_int_sides, 'n_rules':n_rules, 'dim':dim}

    Returns:
        - tuple: origin of tile2 to-be-attached with tile1

    Raises:
        - ValueError: Invalid index or dimension.
    """
    dim = params['dim']
    if dim < 3 and index > 3: raise ValueError("Invalid index. Expected index between 0 and 3 for 2D")

    offsets = [(0, 1, 0), (1, 0, 0), (0, -1, 0), (-1, 0, 0), (0, 0, -1), (0, 0, 1)]
    try:
        return tuple(o + d for o, d in zip(t1_origin, offsets[index])) # origin of tile2 wrt tile1
    except IndexError:
        raise ValueError("Invalid index. Expected index between 0 and 5")
    
#==========================================================================================================#
def zero_sides(t2_origin, tile_coord, rot_interfaces, available_sides):
    """  
    After tile-2 is attached to a tile-1, zero the touching faces of indirectly-blocked neighbouring tiles.

    Args:
        - t2_origin (tuple): origin of tile-2.
        - tile_coord (list of tuples): list of origins of all placed tiles (tile-1(s)).
        - rot_interfaces (list): rotated(if rqrd) interfaces of tile-2.
        - available_sides (list): list of all "open" sides eligible for attachment.

    Returns:
        - rot_interfaces (list): zeroed interfaces of tile2 which are blocked.
        - available_sides (list): zeroed "open" sides which are now blocked.
    """
    (x,y,z) = t2_origin
    if (x,y+1,z) in tile_coord: # if Above tile is blocked
        rot_interfaces[0] = '00' # North edge of tile-2
        tile_num = tile_coord.index((x,y+1,z)) 
        available_sides[6*tile_num+2] = '00' # South edge of above tile
    if (x+1,y,z) in tile_coord: # if Right tile is blocked
        rot_interfaces[1] = '00' # East edge of tile-2
        tile_num = tile_coord.index((x+1,y,z))
        available_sides[6*tile_num+3] = '00' # West edge of right tile
    if (x,y-1,z) in tile_coord:  # if Below tile is blocked
        rot_interfaces[2] = '00' # South edge of tile-2
        tile_num = tile_coord.index((x,y-1,z))
        available_sides[6*tile_num+0] = '00' # North edge of below tile
    if (x-1,y,z) in tile_coord: # if Left tile is blocked
        rot_interfaces[3] = '00' # West edge of tile-2
        tile_num = tile_coord.index((x-1,y,z))
        available_sides[6*tile_num+1] = '00' # East edge of left tile
    if (x,y,z-1) in tile_coord: # if behind tile is blocked
        rot_interfaces[4] = '00' # back edge of tile-2
        tile_num = tile_coord.index((x,y,z-1))
        available_sides[6*tile_num+5] = '00' # front side of behind tile
    if (x,y,z+1) in tile_coord: # if front tile is blocked
        rot_interfaces[5] = '00' # front edge of tile-2
        tile_num = tile_coord.index((x,y,z+1))
        available_sides[6*tile_num+4] = '00' # back side of front tile
    
    return rot_interfaces, available_sides

#==========================================================================================================#
def path_finder(available_sides, available_orient, all_sides, relevant_indices, parent, lineage, choice_line, choice_tree, t2index, comp_tile2_indices, rules_dict):
    """ 
    Generates a choice tree while tracing a particular path (lineage).

    Args:
        - available_sides (list): List of all placed/open sides (with '00' representing 'closed' sides).
        - available_orient (list): List of all orientations corresponding to open sides.
        - all_sides (list): Pool of all tile-2 sides.
        - relevant_indices (list): Indices of "open" sides of placed tiles (tile-1).
        - parent (list): Parent choice line/path, which needs to be further explored now.
        - lineage (list): Stores the current choice line/path.
        - choice_line (list): Path to follow (god's will).
        - choice_tree (list): A tree of all turning points and choices made in this iteration.
        - t2index (int): Index of complementary side chosen from possible complementary tile 2 indices.
        - comp_tile2_indices (list): Indices in all_sides where complementary side is present.
        - rules_dict (dict): Rules of interface attachment (1 <-> 2, 3 <-> 4 etc)

    Returns:
        - random_index (int): Index of the randomly chosen open side of tile-1.
        - t1_orient (str): Orientation of the chosen random side of tile-1.
        - t2index (int): Chosen index of tile-2 from complementary tile 2 indices.
        - lineage (list): Updated lineage with choices made in this iteration.
        - choice_line (list): Updated choice line to be followed.
        - choice_tree (list): Updated choice tree with all possible paths and choices.
        """
    #p#print("available_sides {} \nall_sides {} \nrelevant_indices {} \nparent {} \nlineage {} \nchoice_line {} \nchoice_tree {} \nt2index {} \ncomp_tile2_indices {}".format(available_sides, all_sides, relevant_indices, parent, lineage, choice_line, choice_tree,t2index, comp_tile2_indices))

    # Import packages
    import random
    from utils import extract_underscore

    #?======================= CHOOSE OPEN SIDE (random_index) AND COMPLEMENTARY INDEX (t2index) =======================#
    # -------------------- Tile1 If choice_line is passed, use that to pick open_side else pick randomly --------------#
    if len(choice_line) == 0: random_index = random.sample(relevant_indices,1)[0] #relevant_indices[0] 
    else: random_index = extract_underscore(choice_line[0],'first')
    random_side = available_sides[random_index] # return side at random index (tile1)
    t1_orient = available_orient[random_index] # return orientation of the chosen random side
    other_tile1_indices = list(filter(lambda ind: (ind != random_index), relevant_indices))
    #p#print("Random index = {}, Random side = {}".format(random_index, random_side)) 
    #p#print_msg_box("random index and orientation = "+str(random_index)+" "+str(t1_orient))

    # ------------------ Tile2  If choice_line is passed, use that to pick t2index else pick randomly -----------------#
    comp_side = rules_dict[random_side] # return complementary side of the random side
    comp_tile2_indices = list(filter(lambda pos:(all_sides[pos]==comp_side), range(0, len(all_sides)))) #allcompindices

    if len(choice_line)==0 or choice_line[0][-1] == "_": # If choice line doesnt tell choice of t2index 
        t2index = random.sample(comp_tile2_indices,1)[0] #comp_tile2_indices[0] 
    else: 
        t2index = extract_underscore(choice_line[0], 'second') 
    if len(choice_line) > 0: del choice_line[0]  
    other_tile2_indices = list(filter(lambda ind: (ind != t2index), comp_tile2_indices))
    #p#print("Complementary side = {} \nComplementary tile2 indices = {} \nChosen t2 index {}".format(comp_side, comp_tile2_indices, t2index))
    #p#print_msg_box("t2 index = "+str(t2index))

    lineage.append(str(random_index)+"_"+str(t2index)) # Store choices made in this life
    
    #?================================================= CHOICE TREE ===================================================#
    choice_tree.remove(parent) # Remove parent (only to add it again at the end of list)
    choice_tree.append(parent) # Add parents (at the end) 
    #p#print("Update lineage {} \n Choice_tree with replication {}".format(lineage, choice_tree))

    # Now add all alternate choices to divided parent choice eg if 2 splits into 3,4; write [2,3],[2,4] and add t2index
    temp_choice_tree = []
    for item in choice_tree:
        if item == parent:
            # alternate choices of open side and tile-2 index
            for tile1_index in other_tile1_indices: # as many t1 open side choices
                temp_choice_tree.append(item+[str(tile1_index)+"_"])
            for alt_t2 in other_tile2_indices: # as many t2 indices are possible
                temp_choice_tree.append(item+[str(random_index)+"_"+str(alt_t2)])
            # choices made in this run/life
            temp_choice_tree.append(item+[str(random_index)+"_"+str(t2index)])
        else: temp_choice_tree.append(item) # For paths that are untouched in this life, leave them as they are

    choice_tree = temp_choice_tree.copy()
    #p#print("Temp Choice Tree {} \nChoice Tree {}".format(temp_choice_tree, choice_tree))
    return random_index, t1_orient, t2index, lineage, choice_line, choice_tree

