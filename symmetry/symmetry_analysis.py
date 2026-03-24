""" 
Contains symmetry analysis/utility functions for polycubes.

Functions:
    - permitted_operations(coord_list): Calculate dimensions of the bounding box and return the permitted rotations, reflections and rotoreflections.
    - get_symmetry_order(coord_list): Given a polycube, return the order of its symmetry group after permorming all "permitted" operations.
    - save_symmetry_orders(path, filepath): Given a shape file, calculate the symmetry order of each shape and save it in a text file.
    - get_symmetry_class(sym_operations): Return the class of the polycube (out of 33 classes of O3 group) given the symmetry operations.
    - get_lunnon_data(sorted=True): Return Lunnon's symmetry classes with corresponding orders.
    - construct_3Dclass_and_order_matrices(path, filepath): Construct Lunnon's symmetry class matrix and symmetry order vs size matrix for randomly sampled data. 

"""
# Import packages
import math
import numpy as np
from utils import get_bounding_box, shift_coordinates, convert_tilecoord_to_2d_matrix
from symmetry import reflect_across_plane, rodrigues_rotation, inversion

def permitted_operations(coord_list):
    """ 
    Calculate dimensions of the bounding box and return the permitted rotations, reflections and rotoreflections.

    [Author: Prarthana Agrawal]
    [Date: 12 Feb 2025]
    [Version: 0.1]

    Args:
        - coord_list (list): list of tile-coordinates of the polycube.

    Returns:
        - rot_axes (list): list of axes of rotation.
        - rot_angles (list): list of angles of rotation.
        - ref_planes (list): list of normal vectors of reflection planes.
        - d_list (list): list of distances of reflection planes from the origin.
        - rotoJ_axes (list): list of axes of rotation for rotoreflection.
        - rotoJ_angles (list): list of angles of rotation for rotoreflection.
        - rotoJ_planes (list): list of normal vectors of reflection planes for rotoreflection.
        - rotoJ_d_list (list): list of distances of reflection planes from the origin for rotoreflection.
        - operation_num (list): list of operation numbers (1 to 48).

    Description:
        Not all 48 operations apply on a bounding box since the bbox can be a cube or a cuboid. 
        #! Important update: after debugging, the bbox was forced to be a cube so that all 48 operations apply.
        Only those operations are "permitted" which maps a polycube to itself on a 3D grid. 
        For example, a cuboid may not map to a cuboid under a certain operation.
        - Cube: Octahedral symmetry (48)
        - Rectangular cuboid: D4h Tetragonal symmetry (16) #!CHECK 
        - Irregular cuboid: D2h (8) #! CHECK

    Note: Please refer to the list of symmetry operations (1 to 48) in the notes for more details.
    """
    pi = math.pi # pi value

    # shift coord_list to non-negative values
    coord_list = shift_coordinates(coord_list)

    #? Inversion and Identity are common to all cases, so they are omitted here        

    # Dimensions of bounding box
    bbox_array = np.array(get_bounding_box(coord_list))
    (xmin, ymin, zmin) = np.min(bbox_array, axis=0)
    (xmax, ymax, zmax) = np.max(bbox_array, axis=0)
    #print("Bbox dim =", xmax,ymax,zmax)

    # Midpoints along each axis
    x_mp = (xmin + xmax)/2; y_mp = (ymin + ymax)/2; z_mp = (zmin + zmax)/2
    
    # Length of face diagonals
    norm_yz = math.sqrt(ymax**2 + zmax**2)
    norm_xz = math.sqrt(xmax**2 + zmax**2)
    norm_xy = math.sqrt(xmax**2 + ymax**2)
    
    #? =================================== Irregular cuboid (all sides unequal) ======================================#
    if (xmax != ymax and ymax != zmax and zmax != xmax): # D2h order 8

        rot_axes = list([(1,0,0),(0,1,0),(0,0,1)]) # Face axes: X, Y, Z
        rot_angles = list([pi, pi, pi])
        ref_planes = list([(0,0,1),(1,0,0),(0,1,0)]) # Face planes: XY, YZ, XZ
        d_list = list([z_mp, x_mp, y_mp]) # distance of plane from the origin (face planes pass through midpoints)
        # Operation J = A x E 
        rotoJ_axes = []; rotoJ_angles = [] # roto refers to rotoreflection
        rotoJ_planes = []; rotoJ_d_list = []
        operation_num = [7,8,9,24,25,26,47]

    #? ================================== Rectangular prism (two sides equal) =======================================#
    elif (xmax == ymax and xmax != zmax) or (xmax == zmax and xmax != ymax) or (ymax == zmax and ymax != xmax):
        
        rot_angles = list([
                pi/2, 3*(pi/2),
                pi, pi, pi,
                pi, pi]) # face axis (perpendicular 90,270) // face axis (X,Y,Z 180) // Edge axis (180)
        
        #* ------------------------- X = Y only ------------------------ #
        if (xmax == ymax and xmax != zmax):
            rot_axes = list([
                (0,0,1),(0,0,1),
                (1,0,0),(0,1,0),(0,0,1),
                (xmax,ymax,0), (-xmax,ymax,0)]) # Z, Z // X,Y,Z // OC-DG, AF-BE
            ref_planes = list([
                (0,0,1),(1,0,0),(0,1,0),
                (-xmax, ymax, 0), (xmax, ymax, 0)]) # XY, YZ, XZ // OC-DG, AF-BE
            d_list = list([
                z_mp, x_mp, y_mp,
                0, norm_xy/2])
            # J = A x E
            rotoJ_axes = [(0,0,1),(0,0,1)]
            rotoJ_angles = [pi/2, 3*(pi/2)]
            rotoJ_planes = list([(0,0,1),(0,0,1)]) # XY, XY
            rotoJ_d_list = list([z_mp, z_mp])
            operation_num = [5,6,7,8,9,12,14,24,25,26,29,32,45,46,47] 

        #* ------------------------ X = Z only ----------------------- #
        elif (xmax == zmax and xmax != ymax):
            rot_axes = list([
                (0,1,0),(0,1,0),
                (1,0,0),(0,1,0),(0,0,1),
                (xmax,0,zmax), (-xmax,0,zmax)]) # Y, Y, X, Y, Z, OB-FG, AD-CE
            ref_planes = list([
                (0,0,1),(1,0,0),(0,1,0),
                (-xmax, 0, zmax),(xmax, 0, zmax)]) # XY, YZ, XZ // OB-FG, AD-CE
            d_list = list([
                z_mp, x_mp, y_mp,
                0, norm_xz/2])    
            # Operation J = A x E 
            rotoJ_axes = [(0,1,0),(0,1,0)]
            rotoJ_angles = [pi/2, 3*(pi/2)]
            rotoJ_planes = list([(0,1,0),(0,1,0)]) # XZ, XZ
            rotoJ_d_list = list([y_mp, y_mp])
            operation_num = [3,4,7,8,9,11,13,24,25,26,28,31,43,44,47] 

        #* ------------------------ Y = Z only ---------------------- #
        elif (ymax == zmax and ymax != xmax):
            rot_axes = list([
                (1,0,0),(1,0,0),
                (1,0,0),(0,1,0),(0,0,1),
                (0,ymax,zmax), (0,-ymax,zmax)]) # X, X, X, Y, Z, OA-EG, BD-CF
            ref_planes = list([
                (0,0,1),(1,0,0),(0,1,0),
                (0, -ymax, zmax), (0, ymax, zmax)]) # XY, YZ, XZ // OA-EG, BD-CF
            d_list = list([
                z_mp, x_mp, y_mp,
                0, norm_yz/2])
            # Operation J = A x E 
            rotoJ_axes = [(1,0,0),(1,0,0)]
            rotoJ_angles = [pi/2, 3*(pi/2)]
            rotoJ_planes = list([(1,0,0),(1,0,0)]) # YZ, YZ
            rotoJ_d_list = list([x_mp, x_mp])
            operation_num = [1,2,7,8,9,10,15,24,25,26,27,30,41,42,47]

        else:
            raise Exception("It is not a rectangular prism. Please check dimensions of bbox.")
        
    #? ================================= Cube (all sides equal) ==================================================
    else:
        operation_num = list(np.arange(1,48))

        rot_axes = list([
        (1,0,0),(1,0,0),(0,1,0),(0,1,0),(0,0,1),(0,0,1),(1,0,0),(0,1,0),(0,0,1),
        (0,ymax,zmax),(xmax,0,zmax),(xmax,ymax,0),(-xmax,0,zmax),(-xmax,ymax,0),(0,-ymax,zmax), 
        (xmax, ymax, zmax), (xmax, ymax, zmax), (-xmax, ymax, zmax), (-xmax, ymax, zmax), (xmax, -ymax, zmax), (xmax, -ymax, zmax), (xmax, ymax, -zmax), (xmax, ymax, -zmax)]) 
        # Face axes: X, X, Y, Y, Z, Z, X, Y, Z
        # Edge axes: OA-EG, OB-FG, OC-DG, AD-CE, AF-BE, BD-CF
        # Vertex axes: OG, OG, AE, AE, BF, BF, CD, CD

        rot_angles = list([
        pi/2, 3*(pi/2), pi/2, 3*(pi/2), pi/2, 3*(pi/2), pi, pi, pi, 
        pi, pi, pi, pi, pi, pi, 
        2*(pi/3), 4*(pi/3), 2*(pi/3), 4*(pi/3), 2*(pi/3), 4*(pi/3), 2*(pi/3), 4*(pi/3)])

        ref_planes = list([
        (0,0,1),(1,0,0),(0,1,0),
        (0, -ymax, zmax), (-xmax, 0, zmax), (-xmax, ymax, 0), (0, ymax, zmax), (xmax, 0, zmax), (xmax, ymax, 0)])
        # face planes XY, YZ, XZ
        # edge planes OA-EG, OB-FG, OC-DG, BD-CF, AD-CE, AF-BE

        d_list = list([
            z_mp, x_mp, y_mp,
            0, 0, 0, norm_yz/2, norm_xz/2, norm_xy/2
        ])

       #? Operation J = A x E =========================================================================================
        rotoJ_axes = [(1,0,0),(1,0,0),(0,1,0),(0,1,0),(0,0,1),(0,0,1)]
        rotoJ_angles = [pi/2, 3*(pi/2), pi/2, 3*(pi/2), pi/2, 3*(pi/2)]
        rotoJ_planes = list([(1,0,0),(1,0,0),(0,1,0),(0,1,0),(0,0,1),(0,0,1)]) # YZ, YZ, XZ, XZ, XY, XY
        rotoJ_d_list = list([x_mp, x_mp, y_mp, y_mp, z_mp, z_mp])

    return rot_axes, rot_angles, ref_planes, d_list, rotoJ_axes, rotoJ_angles, rotoJ_planes, rotoJ_d_list, operation_num

#======================================================================================================================#
def get_symmetry_order(coord_list, dim, return_sym_class=False):
    """
    Given a polycube, return the order of its symmetry group after permorming all "permitted" operations.
    
    [Author: Prarthana Agrawal]
    [Date: 12 Feb 2025]
    [Version: 0.1]

    Args:
        - coord_list (list): list of tile-coordinates of the polycube.
        - return_sym_class (bool): return the class of the polycube (out of 33 classes of O3 group)
        - dim (int): Dimensionality of the polycube (2D or 3D). Default is 3D. 
        For 2D, first convert to 2D matrix using convert_tilecoord_to_2d_matrix() and then calculated the symmetry order={1,2,4,8} and class={C1,C2,C4,D1,D2,D4}.

    Returns:
        - sym_count (int): order of the symmetry group of the polycube.

    Raises:
        - None
    """

    # shift coord_list to non-negative values
    coord_list = shift_coordinates(coord_list)

    if dim == 3:

        # Get permitted operations for the polycube's bounding box
        rot_axes, rot_angles, ref_planes, d_list,\
            rotoJ_axes, rotoJ_angles, rotoJ_planes, rotoJ_d_list, operation_num = permitted_operations(coord_list)
        coord_list = [tuple(float(x) for x in tup) for tup in coord_list]

        sym_count = 1 # count number of symmetries
        i = 1 # operation number (matched with notebook polyominoes#3)
        sym_operations = []

        #? ============================================== Rotations ========================================================
        #print("============== Rotations ===============")
        for axis, theta in zip(rot_axes, rot_angles):
            k = axis / np.linalg.norm(axis) # axis of rotation
            rot_vertices = list([rodrigues_rotation(np.array((x, y, z)), k, theta) for (x, y, z) in coord_list])
            shifted_coord = shift_coordinates(rot_vertices)
            if sorted(coord_list) == sorted(shifted_coord): 
                sym_count+=1
                sym_operations.append(i)
            i+=1

        #? ============================================= Reflections =======================================================
        #print("============== Reflections =============")
        for normal, d in zip(ref_planes, d_list):
            normal = normal / np.linalg.norm(normal) # normal of plane
            ref_vertices = list([reflect_across_plane(np.array((x, y, z)), normal, d) for (x, y, z) in coord_list])
            shifted_coord = shift_coordinates(ref_vertices)
            shifted_coord = list([tuple(np.around((x,y,z), 4)) for (x,y,z) in shifted_coord])
            if sorted(coord_list) == sorted(shifted_coord): 
                sym_count+=1
                sym_operations.append(i)
            i+=1

        #? ============================================ Rotoreflections ====================================================
        if len(rot_axes) == 23: # this means the bounding box is a cube
            #print("============== H = D x K ===========")
            # Operation H = D x K
            vertex_axes = rot_axes[15:]
            vertex_angles = rot_angles[15:]
            for axis, theta in zip(vertex_axes, vertex_angles):
                k = axis / np.linalg.norm(axis) # axis of rotation
                rot_vertices = list([rodrigues_rotation(np.array((x, y, z)), k, theta) for (x, y, z) in coord_list])
                inverted_coord = inversion(rot_vertices)
                shifted_coord = shift_coordinates(inverted_coord)
                if sorted(coord_list) == sorted(shifted_coord): 
                    sym_count+=1
                    sym_operations.append(i)
                i+=1

        #? Operation J = A x E =========================================================================================
        if len(rot_axes) != 3: # this means the bounding box is either a cube or a rectangular prism (two sides equal)
            #print("============== J = A x E ===========")
            for j in range(len(rotoJ_axes)):
                # rotation A
                axis = rotoJ_axes[j]; theta = rotoJ_angles[j]
                k = axis / np.linalg.norm(axis) # axis of rotation
                rot_vertices = list([rodrigues_rotation(np.array((x, y, z)), k, theta) for (x, y, z) in coord_list])
                
                # roto reflection A x E
                normal = rotoJ_planes[j]; d = rotoJ_d_list[j]
                normal = normal / np.linalg.norm(normal) # normal of plane #! Not really necessary here
                ref_vertices = list([reflect_across_plane(np.array((x, y, z)), normal, d) for (x, y, z) in rot_vertices])
                shifted_coord = shift_coordinates(ref_vertices)
                shifted_coord = list([tuple(np.around((x,y,z), 4)) for (x,y,z) in shifted_coord])
                if sorted(coord_list) == sorted(shifted_coord): 
                    sym_count+=1
                    sym_operations.append(i)
                i+=1

        #? Inversion ====================================================================================================
        ref_vertices = inversion(coord_list)
        shifted_coord = shift_coordinates(ref_vertices)
        shifted_coord = list([tuple(np.around((x,y,z), 4)) for (x,y,z) in shifted_coord])
        if sorted(coord_list) == sorted(shifted_coord): 
            sym_count+=1
            sym_operations.append(i)
        i+=1
            
        sym_operations.append(i) # identity operation i= 48 
        #print("Order =", sym_count)

        if return_sym_class:
            sym_class = get_symmetry_class(sym_operations)
        else:
            sym_class = None

    elif dim == 2:
        sym_count, sym_class = get_2d_symmetry_order_and_class(coord_list)

    else: 
        raise Exception("Dimension not supported. Please use dim=2 or dim=3.")
    #print("Symmetry operations:", sym_operations)
    return sym_count, sym_class

#======================================================================================================#
def get_2d_symmetry_order_and_class(coords):
    """
    Given a polyomino, return the order of its symmetry group after permorming all square symmetry operations.
    """ 
    # Convert tile coordinates to 2D matrix
    shape_matrix = convert_tilecoord_to_2d_matrix(coords)

    # 8 symmetry operations in 2D
    A = shape_matrix
    ACW90 = np.rot90(A)
    ACW180 = np.rot90(ACW90)
    ACW270 = np.rot90(ACW180)
    R_lr = np.fliplr(A)
    R_ud = np.flipud(A)
    R_d1 = np.transpose(A)
    R_d2 = np.flipud(np.fliplr(np.transpose(A)))

    operations = [ACW90, ACW180, ACW270, R_lr, R_ud, R_d1, R_d2]
    invariant = list([1])
    count_sym_order = 1

    for i, op in enumerate(operations):
        
        if np.array_equal(op, A) == True: 
            invariant.append(i+2)
            count_sym_order += 1

    if invariant == [1]: sym = 'C1'
    elif invariant == [1,3]: sym = 'C2'
    elif invariant == [1,2,3,4]: sym = 'C4'
    elif invariant == [1,5] or invariant == [1,6] or invariant == [1,7] or invariant == [1,8]: sym = 'D1'
    elif invariant == [1,3,5,6] or invariant == [1,3,5,7] or invariant == [1,3,5,8] or invariant == [1,3,6,7] or invariant == [1,3,6,8] or invariant == [1,3,7,8]: sym = 'D2' 
    elif invariant == [1,2,3,4,5,6,7,8]: sym = 'D4'
    else: raise Exception("Symmetry Unknown")

    return count_sym_order, sym


#? ===================================================================================================================#

def save_symmetry_orders(path, filepath, dim, save_sym_class = False,):
    """
    Given a shape file, calculate the symmetry order of each shape and save it in a text file.
    
    Args:
        - path (str): path to the directory containing the shape file.
        - filepath (str): path of the shape file.
        - save_sym_class (bool): save the symmetry class of the polycube (out of 33 classes of O3 group)

    Returns:
        - sym_orders (list): list of symmetry orders.
        - sym_order_file (str): text file containing symmetry orders.
        - sym_class_file (str): text file containing symmetry classes.
    """
    if dim == 3:
        
        shape_file = open(path + filepath + 'valid_shapes.txt', 'r')
        sym_orders = [] # list of symmetry orders
        sym_classes = []
        for line in shape_file:
            coords = eval(line.strip())
            shift_coords = shift_coordinates(coords) # not strictly necessary, already implemented in get_symmetry_order
            sym_order, sym_class = get_symmetry_order(shift_coords, dim, return_sym_class=save_sym_class)
            sym_orders.append(sym_order)
            sym_classes.append(sym_class)

    elif dim == 2:

        shape_file = open(path + filepath + 'valid_shapes.txt', 'r')
        sym_orders = [] # list of symmetry orders
        sym_classes = []
        for line in shape_file:
            coords = eval(line.strip())
            sym_order, sym_class = get_2d_symmetry_order_and_class(coords)
            sym_orders.append(sym_order)
            sym_classes.append(sym_class)

    else: raise Exception("Dimension not supported. Please use dim=2 or dim=3.")
        
    # Save symmetry orders in a text file
    sym_order_file = open(path + filepath + 'symmetry_orders.txt', 'w')
    for order in sym_orders:
        sym_order_file.write(str(order) + '\n')
    sym_order_file.close()

    if save_sym_class:
        # Save symmetry classes in a text file
        sym_class_file = open(path + filepath + 'symmetry_classes.txt', 'w')
        for sym_class in sym_classes:
            sym_class_file.write(str(sym_class) + '\n')
        sym_class_file.close()

#======================================================================================================================#

def get_symmetry_class(sym_operations):
    """ 
    Return the class of the polycube (out of 33 classes of O3 group) given the symmetry operations.

    Args:
        - sym_operations (list): list of symmetry operations (from 1 to 48, please note that inversion is 47 and identity is 48).
    
    Returns:
        - sym_class ('str'): class of the polycube.

    References:  W.F.Lunnon: Symmetry of cubical and general polyominoes, pp. 101-108 of R. C. Read, editor, Graph Theory and Computing. Academic Press, NY, 1972.
    """
    # mapping of symmetry operations to symmetry classes (Lunnon)
    symmetry_classes_dict = {'A':[1,2,3,4,5,6], 'B':[7,8,9], 'C':[10,11,12,13,14,15], 'D':[16,17,18,19,20,21,22,23], 'E':[24,25,26], 'F':[27,28,29,30,31,32], 'H':[33,34,35,36,37,38,39,40],'J':[41,42,43,44,45,46],'K':[47],'I':[48]}

    # order of the symmetry group
    order = len(sym_operations)

    if order == 1:
        return 'I'
    
    if order == 2:
        for key in ['B', 'C', 'E', 'F', 'K']:
            if len(set(sym_operations) & set(symmetry_classes_dict[key])) == 1:
                return key

    if order == 3:
        if len(set(sym_operations) & set(symmetry_classes_dict['D'])) == 2:
            return 'D'
    
    if order == 4:
        conditions = {
            'A': ('B', 1, 'A', 2),
            'J': ('B', 1, 'J', 2),
            'BB': ('B', 3),
            'BC': ('B', 1, 'C', 2),
            'BE': ('B', 1, 'E', 1, 'K', 1),
            'BF': ('B', 1, 'F', 2),
            'CE': ('C', 1, 'E', 1, 'F', 1),
            'CK': ('C', 1, 'F', 1, 'K', 1),
            'EE': ('B', 1, 'E', 2)
        }
        
        for sym_class, cond in conditions.items():
            if all(len(set(sym_operations) & set(symmetry_classes_dict[c])) == n for c, n in zip(cond[::2], cond[1::2])):
                return sym_class
    
    if order == 6:
        conditions = {
            'H': ('K', 1, 'D', 2, 'H', 2),
            'CD': ('D', 2, 'C', 3),
            'FF': ('D', 2, 'F', 3)
        }

        for sym_class, cond in conditions.items():
            if all(len(set(sym_operations) & set(symmetry_classes_dict[c])) == n for c, n in zip(cond[::2], cond[1::2])):
                return sym_class
    
    if order == 8:
        conditions = {
            'AB': ('A', 2, 'C', 2, 'B', 3),
            'AE': ('B', 1, 'E', 1, 'K', 1, 'A', 2, 'J', 2),
            'BFF': ('F', 2, 'J', 2, 'B', 3),
            'CJ': ('B', 1, 'C', 2, 'E', 2, 'J', 2),
            'EEE': ('K', 1, 'B', 3, 'E', 3),
            'EF': ('B', 1, 'A', 2, 'E', 2, 'F', 2),
            'EFF': ('B', 1, 'E', 1, 'K', 1, 'C', 2, 'F', 2)
        }
        
        for sym_class, cond in conditions.items():
            if all(len(set(sym_operations) & set(symmetry_classes_dict[c])) == n for c, n in zip(cond[::2], cond[1::2])):
                return sym_class
    
    if order == 12:
        conditions = {
            'BD': ('B', 3, 'D', 8),
            'CF': ('K', 1, 'D', 2, 'H', 2, 'C', 3, 'F', 3)
        }
        
        for sym_class, cond in conditions.items():
            if all(len(set(sym_operations) & set(symmetry_classes_dict[c])) == n for c, n in zip(cond[::2], cond[1::2])):
                return sym_class
    
    if order == 16:
        conditions = {
            'BBC': ('K', 1, 'A', 2,'C', 2, 'F', 2, 'J', 2,'B', 3, 'E', 3)
        }
        for sym_class, cond in conditions.items():
            if all(len(set(sym_operations) & set(symmetry_classes_dict[c])) == n for c, n in zip(cond[::2], cond[1::2])):
                return sym_class
    
    if order == 24:
        conditions = {
            'CCC': ('B', 3, 'F', 6, 'J', 6, 'D', 8),
            'DEE': ('K', 1, 'B', 3, 'E', 3, 'D', 8, 'H', 8),
            'R': ('B', 3, 'A', 6, 'C', 6, 'D', 8)
        }
        
        for sym_class, cond in conditions.items():
            if all(len(set(sym_operations) & set(symmetry_classes_dict[c])) == n for c, n in zip(cond[::2], cond[1::2])):
                return sym_class
    
    if order == 48:
        return 'G'
    
    return None
     

#======================================================================================================================#

def get_lunnon_data(sort=True):
    """ 
    Retrieve all 33 symmetry classes with their orders and class structures.

    Args:
        - sort (bool): If True, sort the symmetry classes, primarily in increasing symmetry order. Incase of ties, sort in increasing unique number of symmetry operations (class structure with no duplicates).
    
    Returns:
        - lunnon_sym_classes_order (dict): dictionary keys: symmetry classes, values: orders

    Author: Prarthana Agrawal
    Date: 19 Mar 2025
    """
    # Dictionary key: value = class: order
    lunnon_sym_classes_order = dict({'I':1, 'A':4, 'B':2, 'C':2, 'D':3, 'E':2, 'F':2, 'H':6, 'J':4, 'K':2, 'BB':4, 'BC':4,'BE':4,'BF':4,'CE':4,'CK':4,'EE':4,'CD':6,'FF':6,'AB':8,'AE':8,'BFF':8, 'CJ':8,'EEE':8,'EF':8,'EFF':8,'BD':12,'CF':12,'BBC':16,'CCC':24,'DEE':24,'R':24,'G':48})

    if not sort:
        return lunnon_sym_classes_order
    
    # lunnon symmetry class structure (only unique operations, no repetitions eg 2A is written as A)
    # Dictionary key: value = class: operations
    lunnon_sym_class_structure = {
        'I': ['I'], 'A': ['I', 'B', 'A'], 'B': ['I', 'B'], 'C': ['I', 'C'], 'D': ['I', 'D'], 'E': ['I', 'E'], 'F': ['I', 'F'],
        'H': ['I', 'K', 'D', 'H'], 'J': ['I', 'B', 'J'], 'K': ['I', 'K'], 'BB': ['I', 'B'], 'BC': ['I', 'B', 'C'],
        'BE': ['I', 'B', 'E', 'K'], 'BF': ['I', 'B', 'F'], 'CE': ['I', 'C', 'E', 'F'], 'CK': ['I', 'C', 'F', 'K'], 
        'EE': ['I', 'B', 'E'], 'CD': ['I', 'D', 'C'], 'FF': ['I', 'D', 'F'],
        'AB': ['I', 'A', 'C', 'B'], 'AE': ['I', 'B', 'E', 'K', 'A', 'J'], 
        'BFF': ['I', 'F', 'J', 'B'], 'CJ': ['I', 'B', 'C', 'E', 'J'], 
        'EEE': ['I', 'K', 'B', 'E'], 'EF': ['I', 'B', 'A', 'E', 'F'], 
        'EFF': ['I', 'B', 'E', 'K', 'C', 'F'], 'BD': ['I', 'B', 'D'], 
        'CF': ['I', 'K', 'D', 'H', 'C', 'F'], 
        'BBC': ['I', 'K', 'A', 'C', 'F', 'J', 'B', 'E'], 
        'CCC': ['I', 'B', 'F', 'J', 'D'], 
        'DEE': ['I', 'K', 'B', 'E', 'D', 'H'], 
        'R': ['I', 'B', 'A', 'C', 'D'], 
        'G': ['I', 'K', 'B', 'E', 'A', 'C', 'F', 'J', 'D', 'H']
    }

    #sorted_lunnon_sym_classes_order has elements ordered first by increasing symmetry order and then, for ties, by increasing the number of operations.
    sorted_lunnon_sym_classes_order = dict(sorted(
        lunnon_sym_classes_order.items(),
        key=lambda item: (item[1], len(lunnon_sym_class_structure[item[0]]))
    ))

    lunnon_sym_classes_order = sorted_lunnon_sym_classes_order

    return lunnon_sym_classes_order

#======================================================================================================================#
def construct_3Dclass_and_order_matrices(path, filepath):
    """ 
    Constructs two matrices using randomly sampled 3D polycube data:
        1. A symmetry class matrix: rows represent Lunnon's symmetry classes, columns represent size `n` (size of polycube = number of cubes).
        2. A symmetry order matrix: rows represent Lunnon's symmetry orders, columns represent size `n`.

    For each valid polycube shape:
        - It determines the shape's size `n`.
        - Retrieves its symmetry class and symmetry order.
        - Increments the respective matrix positions with its frequency.

    Only sizes up to a maximum `n` (e.g. 22) are considered, consistent with OEIS data.

    Args:
        - path (str): path to the directory containing the shape files.
        - filepath (str): relative path of the folder where the shape files are stored.

    Returns:
        - symclass_matrix (np.ndarray): matrix [num_classes x max_size]
        - symorder_matrix (np.ndarray): matrix [num_orders x max_size]
    """
    # Import dependencies
    from collections import defaultdict
    from utils.oeis_data import get_oeis_3Dclass_matrix, get_oeis_3Dorder_matrix

    # Load precomputed data files
    sym_classes_file = np.loadtxt(path + filepath + 'symmetry_classes.txt', dtype='str')
    sym_orders_file = np.loadtxt(path + filepath + 'symmetry_orders.txt')
    frequency_file = np.loadtxt(path + filepath + 'frequency.txt')
    
    # Read valid shapes from file
    with open(path + filepath + 'valid_shapes.txt', 'r') as shape_file:
        shape_content = shape_file.readlines()
    nshapes = len(shape_content)  # Number of unique polycube shapes

    # Get Lunnon's symmetry class→order mapping
    lunnon_classes_orders_dict = get_lunnon_data(sort=True)
    lunnon_sym_classes_list = list(lunnon_classes_orders_dict.keys())

    # Group classes by symmetry order 
    grouped_order_sym_class_dict = defaultdict(list)
    for sym_class, sym_order in lunnon_classes_orders_dict.items():
        grouped_order_sym_class_dict[sym_order].append(sym_class)
    grouped_order_sym_class_dict = dict(sorted(grouped_order_sym_class_dict.items()))
    orders_list = list(grouped_order_sym_class_dict.keys())

    # Get max shape size from OEIS data
    oeis_3Dclass_matrix = get_oeis_3Dclass_matrix()
    max_size = oeis_3Dclass_matrix.shape[1] # max size (n) for which all sym classes data is available on OEIS database

    # Initialize matrices
    symclass_matrix = np.zeros((len(lunnon_sym_classes_list), max_size), dtype=np.int64) # rows: sym classes, cols: size
    symorder_matrix = np.zeros((len(orders_list), max_size), dtype=np.int64) # rows: orders, cols: size

    # Construct both matrices
    """
    # 1. Iterate over each shape
    # 2. Determine size of the shape (number of unit cubes)
    # 3. Retrieve symmetry class and order
    # 4. Identify row indices for class and order matrices
    # 5. Update matrices with frequency of the shape
    """

    for i in range(nshapes):
        shape = eval(shape_content[i])         # convert shape string to list of coordinates
        n = len(shape)                         # size of the shape = number of unit cubes
        
        if n > max_size:
            continue  # skip shapes beyond OEIS-supported size

        freq = frequency_file[i]               # frequency of this shape in dataset
        sym_class = sym_classes_file[i]        # symmetry class of this shape
        sym_order = sym_orders_file[i]         # symmetry order of this shape

        # Identify row indices
        class_row = lunnon_sym_classes_list.index(sym_class)
        order_row = orders_list.index(sym_order)
        
        col = n - 1  # 0-based indexing for shape size

        # Update matrices
        symclass_matrix[class_row, col] += freq
        symorder_matrix[order_row, col] += freq

    return symclass_matrix, symorder_matrix

#======================================================================================================================#
def unique_3Dclass_and_order_matrices(path, filepath):
    """ 
    Constructs two matrices using randomly sampled 3D polycube data (
    #!count only unique valid shapes not their frequencies)
        1. A symmetry class matrix: rows represent Lunnon's symmetry classes, columns represent size `n` (size of polycube = number of cubes).
        2. A symmetry order matrix: rows represent Lunnon's symmetry orders, columns represent size `n`.

    For each valid polycube shape:
        - It determines the shape's size `n`.
        - Retrieves its symmetry class and symmetry order.
        - Increments the respective matrix positions by 1.

    Only sizes up to a maximum `n` (e.g. 22) are considered, consistent with OEIS data.

    Args:
        - path (str): path to the directory containing the shape files.
        - filepath (str): relative path of the folder where the shape files are stored.

    Returns:
        - symclass_matrix (np.ndarray): matrix [num_classes x max_size]
        - symorder_matrix (np.ndarray): matrix [num_orders x max_size]
    """
    # Import dependencies
    from collections import defaultdict
    from utils.oeis_data import get_oeis_3Dclass_matrix, get_oeis_3Dorder_matrix

    # Load precomputed data files
    sym_classes_file = np.loadtxt(path + filepath + 'symmetry_classes.txt', dtype='str')
    sym_orders_file = np.loadtxt(path + filepath + 'symmetry_orders.txt')
    frequency_file = np.loadtxt(path + filepath + 'frequency.txt')
    
    # Read valid shapes from file
    with open(path + filepath + 'valid_shapes.txt', 'r') as shape_file:
        shape_content = shape_file.readlines()
    nshapes = len(shape_content)  # Number of unique polycube shapes

    # Get Lunnon's symmetry class→order mapping
    lunnon_classes_orders_dict = get_lunnon_data(sort=True)
    lunnon_sym_classes_list = list(lunnon_classes_orders_dict.keys())

    # Group classes by symmetry order 
    grouped_order_sym_class_dict = defaultdict(list)
    for sym_class, sym_order in lunnon_classes_orders_dict.items():
        grouped_order_sym_class_dict[sym_order].append(sym_class)
    grouped_order_sym_class_dict = dict(sorted(grouped_order_sym_class_dict.items()))
    orders_list = list(grouped_order_sym_class_dict.keys())

    # Get max shape size from OEIS data
    oeis_3Dclass_matrix = get_oeis_3Dclass_matrix()
    max_size = oeis_3Dclass_matrix.shape[1] # max size (n) for which all sym classes data is available on OEIS database

    # Initialize matrices
    symclass_matrix = np.zeros((len(lunnon_sym_classes_list), max_size), dtype=np.int64) # rows: sym classes, cols: size
    symorder_matrix = np.zeros((len(orders_list), max_size), dtype=np.int64) # rows: orders, cols: size

    # Construct both matrices
    """
    # 1. Iterate over each shape
    # 2. Determine size of the shape (number of unit cubes)
    # 3. Retrieve symmetry class and order
    # 4. Identify row indices for class and order matrices
    # 5. Update matrices with frequency of the shape
    """

    for i in range(nshapes):
        shape = eval(shape_content[i])         # convert shape string to list of coordinates
        n = len(shape)                         # size of the shape = number of unit cubes
        
        if n > max_size:
            continue  # skip shapes beyond OEIS-supported size

        freq = frequency_file[i]               # frequency of this shape in dataset
        sym_class = sym_classes_file[i]        # symmetry class of this shape
        sym_order = sym_orders_file[i]         # symmetry order of this shape

        # Identify row indices
        class_row = lunnon_sym_classes_list.index(sym_class)
        order_row = orders_list.index(sym_order)
        
        col = n - 1  # 0-based indexing for shape size

        # Update matrices
        symclass_matrix[class_row, col] += 1 #! Note this step
        symorder_matrix[order_row, col] += 1 #! Note this step

    return symclass_matrix, symorder_matrix


#====================================================================================================#
def sym_classes_per_order(dim):
    """
    According to the dimensionality of polyomino assembly, return the list of symmetry classes corresponding to each symmetry order.

    eg in 3D: {1: ['I'], 2: ['B', 'C', 'E', 'F', 'K'], 3: ['D'], 4: ['BB', 'A', 'J', 'BC', 'BF', 'EE', 'BE', 'CE', 'CK'], 6: ['CD', 'FF', 'H'], 8: ['AB', 'BFF', 'EEE', 'CJ', 'EF', 'AE', 'EFF'], 12: ['BD', 'CF'], 16: ['BBC'], 24: ['CCC', 'R', 'DEE'], 48: ['G']}

    eg in 2D: {1: ['C1'], 2: ['C2', 'D1'], 4: ['C4', 'D2'], 8: ['D4']}
    
    Args:
        - dim (int): dimensionality of the polycube (2D or 3D).

    Returns:
        - sym_classes_per_order_dict (dict): dictionary with keys as symmetry orders and values as list of symmetry classes.
    """
    from collections import defaultdict
    from symmetry import get_lunnon_data

    if dim == 3:

        lunnon_classes_orders_dict = get_lunnon_data(sort=True)
        grouped_order_sym_class_dict = defaultdict(list)
        for key, value in lunnon_classes_orders_dict.items():
            grouped_order_sym_class_dict[value].append(key)
        sym_classes_per_order_dict = dict(sorted(grouped_order_sym_class_dict.items()))

    elif dim == 2:

        sym_classes_per_order_dict = {
            1: ['C1'],
            2: ['C2', 'D1'],
            4: ['C4', 'D2'],
            8: ['D4']
        }

    else:
        raise Exception("Dimension not supported. Please use dim=2 or dim=3.")
    return sym_classes_per_order_dict

#======================================================================================================================#
def get_classes_and_orders_dict(dim):
    """
    Given dim (2D or 3D), return dictionary:
        sym_classes_and_orders_dict: keys as symmetry classes and values as symmetry orders.
    """
    if dim == 3:
        sym_classes_and_orders_dict = dict({'I':1, 'A':4, 'B':2, 'C':2, 'D':3, 'E':2, 'F':2, 'H':6, 'J':4, 'K':2, 'BB':4, 'BC':4,'BE':4,'BF':4,'CE':4,'CK':4,'EE':4,'CD':6,'FF':6,'AB':8,'AE':8,'BFF':8, 'CJ':8,'EEE':8,'EF':8,'EFF':8,'BD':12,'CF':12,'BBC':16,'CCC':24,'DEE':24,'R':24,'G':48})

    elif dim == 2:
        sym_classes_and_orders_dict = dict({'C1':1, 'C2':2, 'D1':2, 'C4':4, 'D2':4, 'D4':8})

    else: 
        raise Exception("Dimension not supported. Please use dim=2 or dim=3.")
    return sym_classes_and_orders_dict


#======================================================================================================================#
def get_symclass_to_vector_dict(dim):
    """ 
    Given a symmetry class and dimensionality, return a vector representing the presence of each symmetry operation in that class. 
    
    """
    symmetry_dict = {
    "I":   {"I":1,"A":0,"B":0,"C":0,"D":0,"E":0,"F":0,"H":0,"J":0,"K":0},
    "A":   {"I":1,"A":2,"B":1,"C":0,"D":0,"E":0,"F":0,"H":0,"J":0,"K":0},
    "B":   {"I":1,"A":0,"B":1,"C":0,"D":0,"E":0,"F":0,"H":0,"J":0,"K":0},
    "C":   {"I":1,"A":0,"B":0,"C":1,"D":0,"E":0,"F":0,"H":0,"J":0,"K":0},
    "D":   {"I":1,"A":0,"B":0,"C":0,"D":2,"E":0,"F":0,"H":0,"J":0,"K":0},
    "E":   {"I":1,"A":0,"B":0,"C":0,"D":0,"E":1,"F":0,"H":0,"J":0,"K":0},
    "F":   {"I":1,"A":0,"B":0,"C":0,"D":0,"E":0,"F":1,"H":0,"J":0,"K":0},
    "H":   {"I":1,"A":0,"B":0,"C":0,"D":2,"E":0,"F":0,"H":2,"J":0,"K":1},
    "J":   {"I":1,"A":0,"B":1,"C":0,"D":0,"E":0,"F":0,"H":0,"J":2,"K":0},
    "K":   {"I":1,"A":0,"B":0,"C":0,"D":0,"E":0,"F":0,"H":0,"J":0,"K":1},
    "BB":  {"I":1,"A":0,"B":3,"C":0,"D":0,"E":0,"F":0,"H":0,"J":0,"K":0},
    "BC":  {"I":1,"A":0,"B":1,"C":2,"D":0,"E":0,"F":0,"H":0,"J":0,"K":0},
    "BE":  {"I":1,"A":0,"B":1,"C":0,"D":0,"E":1,"F":0,"H":0,"J":0,"K":1},
    "BF":  {"I":1,"A":0,"B":1,"C":0,"D":0,"E":0,"F":2,"H":0,"J":0,"K":0},
    "CE":  {"I":1,"A":0,"B":0,"C":1,"D":0,"E":1,"F":1,"H":0,"J":0,"K":0},
    "CK":  {"I":1,"A":0,"B":0,"C":1,"D":0,"E":0,"F":1,"H":0,"J":0,"K":1},
    "EE":  {"I":1,"A":0,"B":1,"C":0,"D":0,"E":2,"F":0,"H":0,"J":0,"K":0},
    "CD":  {"I":1,"A":0,"B":0,"C":3,"D":2,"E":0,"F":0,"H":0,"J":0,"K":0},
    "FF":  {"I":1,"A":0,"B":0,"C":0,"D":2,"E":0,"F":3,"H":0,"J":0,"K":0},
    "AB":  {"I":1,"A":2,"B":3,"C":2,"D":0,"E":0,"F":0,"H":0,"J":0,"K":0},
    "AE":  {"I":1,"A":2,"B":1,"C":0,"D":0,"E":1,"F":0,"H":0,"J":2,"K":1},
    "BFF": {"I":1,"A":0,"B":3,"C":0,"D":0,"E":0,"F":2,"H":0,"J":2,"K":0},
    "CJ":  {"I":1,"A":0,"B":1,"C":2,"D":0,"E":2,"F":0,"H":0,"J":2,"K":0},
    "EEE": {"I":1,"A":0,"B":3,"C":0,"D":0,"E":3,"F":0,"H":0,"J":0,"K":1},
    "EF":  {"I":1,"A":2,"B":1,"C":0,"D":0,"E":2,"F":2,"H":0,"J":0,"K":0},
    "EFF": {"I":1,"A":0,"B":1,"C":2,"D":0,"E":1,"F":2,"H":0,"J":0,"K":1},
    "BD":  {"I":1,"A":0,"B":3,"C":0,"D":8,"E":0,"F":0,"H":0,"J":0,"K":0},
    "CF":  {"I":1,"A":0,"B":0,"C":3,"D":2,"E":0,"F":3,"H":2,"J":0,"K":1},
    "BBC": {"I":1,"A":2,"B":3,"C":2,"D":0,"E":3,"F":2,"H":0,"J":2,"K":1},
    "CCC": {"I":1,"A":0,"B":3,"C":0,"D":8,"E":0,"F":6,"H":0,"J":6,"K":0},
    "DEE": {"I":1,"A":0,"B":3,"C":0,"D":8,"E":3,"F":0,"H":8,"J":0,"K":1},
    "R":   {"I":1,"A":6,"B":3,"C":6,"D":8,"E":0,"F":0,"H":0,"J":0,"K":0},
    "G":   {"I":1,"A":6,"B":3,"C":6,"D":8,"E":3,"F":6,"H":8,"J":6,"K":1},
    }

    class_vectors_dict = dict()
    if dim == 3:
        for sym_class in symmetry_dict:
            sym_vector = list(symmetry_dict[sym_class].values())
            class_vectors_dict[sym_class] = sym_vector
            
    if dim == 2:
        raise Exception("This function is currently implemented for 3D symmetry classes only. Please use dim=3.")
    
    return class_vectors_dict


#======================================================================================================================#
def convert_class_to_vector(sym_class, dim):
    """ 
    Given a symmetry class and dimensionality, return a vector representing the presence of each symmetry operation in that class. 
    """
    symclass_to_vector_dict = get_symclass_to_vector_dict(dim)
    if sym_class not in symclass_to_vector_dict:
        raise Exception(f"Symmetry class {sym_class} not found for dimension {dim}.")
    return symclass_to_vector_dict[sym_class]