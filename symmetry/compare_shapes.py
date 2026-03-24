"""
# Description: Functions to compare two polycubes- whether they are equivalent under any symmetry operation.

# Functions:
    - brute_force_comparison(first_shape, second_shape): performs all 3D operations possible on a polycube.
    - compare_polycubes(tile_coord, valid_shape_i): checks if the new shape is similar to an already encountered shape.

    Difference between the two functions is that compare_polycubes() first checks some elimination criteria before calling brute_force_comparison(). This makes it more efficient than brute_force_comparison().

# Dependencies:
    - numpy
    - math
    - utils.shift_coordinates
    - utils.get_bounding_box
    - utils.return_lengths
    - symmetry.rodrigues_rotation
    - symmetry.reflect_across_plane
    - symmetry.inversion
"""

# Import packages
import math
import numpy as np
from utils import shift_coordinates, get_bounding_box, return_lengths
from symmetry import rodrigues_rotation, reflect_across_plane, inversion

def brute_force_comparison(first_shape, second_shape):
    """ 
    Performs all 3D operations possible on a polycube.
    Check if first_shape matches second_shape under any symmetry operation.
    
    Args:
        - first_shape (list of tuples): list of cube-origins of first polycube.
        - second_shape (list of tuples): list of cube-origins of second polycube.

    Returns:
        - 1 if first_shape is equivalent to second_shape under any symmetry operation.
        - 0 otherwise.
    
    Raises:
        - None
    
    #! TODO: Please create a reference for the operations and their corresponding numbers.

    #! NOTE: Here I am not checking the order of symmetry of a shape. That can be done in post-processing once all shapes are obtained.
    """
    pi = math.pi

    # Make sure all coordinates are non-negative #? THIS IS IMPORTANT!
    first_shape = shift_coordinates(first_shape)
    second_shape = shift_coordinates(second_shape)

    # Identity Operation
    if sorted(first_shape) ==  sorted(second_shape):
        return 1
    
    # Dimensions of bounding box
    bbox_array = np.array(get_bounding_box(second_shape))
    (xmin, ymin, zmin) = np.min(bbox_array, axis=0)
    (xmax, ymax, zmax) = np.max(bbox_array, axis=0)
    #p#print("Bbox dim =", xmax,ymax,zmax)

    # Midpoints along each axis
    x_mp = (xmin + xmax)/2; y_mp = (ymin + ymax)/2; z_mp = (zmin + zmax)/2
    
    # Length of face diagonals
    norm_yz = math.sqrt(ymax**2 + zmax**2)
    norm_xz = math.sqrt(xmax**2 + zmax**2)
    norm_xy = math.sqrt(xmax**2 + ymax**2)
    
    #operation_num = list(np.arange(1,48))

    # Rotational axes
    rot_axes = list([
    (1,0,0),(1,0,0),(0,1,0),(0,1,0),(0,0,1),(0,0,1),(1,0,0),(0,1,0),(0,0,1),
    (0,ymax,zmax),(xmax,0,zmax),(xmax,ymax,0),(-xmax,0,zmax),(-xmax,ymax,0),(0,-ymax,zmax), 
    (xmax, ymax, zmax), (xmax, ymax, zmax), (-xmax, ymax, zmax), (-xmax, ymax, zmax), (xmax, -ymax, zmax), (xmax, -ymax, zmax), (xmax, ymax, -zmax), (xmax, ymax, -zmax)]) 
    # Face axes: X, X, Y, Y, Z, Z, X, Y, Z
    # Edge axes: OA-EG, OB-FG, OC-DG, AD-CE, AF-BE, BD-CF
    # Vertex axes: OG, OG, AE, AE, BF, BF, CD, CD

    # Rotational angles
    rot_angles = list([
    pi/2, 3*(pi/2), pi/2, 3*(pi/2), pi/2, 3*(pi/2), pi, pi, pi, 
    pi, pi, pi, pi, pi, pi, 
    2*(pi/3), 4*(pi/3), 2*(pi/3), 4*(pi/3), 2*(pi/3), 4*(pi/3), 2*(pi/3), 4*(pi/3)])

    # Reflection planes
    ref_planes = list([
    (0,0,1),(1,0,0),(0,1,0),
    (0, -ymax, zmax), (-xmax, 0, zmax), (-xmax, ymax, 0), (0, ymax, zmax), (xmax, 0, zmax), (xmax, ymax, 0)])
    # face planes XY, YZ, XZ
    # edge planes OA-EG, OB-FG, OC-DG, BD-CF, AD-CE, AF-BE

    # Reflection distances
    d_list = list([
        z_mp, x_mp, y_mp,
        0, 0, 0, norm_yz/2, norm_xz/2, norm_xy/2
    ])

    #? Operation J = A x E =========================================================================================
    rotoJ_axes = [(1,0,0),(1,0,0),(0,1,0),(0,1,0),(0,0,1),(0,0,1)]
    rotoJ_angles = [pi/2, 3*(pi/2), pi/2, 3*(pi/2), pi/2, 3*(pi/2)]
    rotoJ_planes = list([(1,0,0),(1,0,0),(0,1,0),(0,1,0),(0,0,1),(0,0,1)]) # YZ, YZ, XZ, XZ, XY, XY
    rotoJ_d_list = list([x_mp, x_mp, y_mp, y_mp, z_mp, z_mp])

    #======================================== Perform operations ======================================================#
    first_shape = [tuple(float(x) for x in tup) for tup in first_shape]
    sorted_first_shape = sorted(first_shape)
    #p#print("First shape", sorted_first_shape)
    
    #i=0 # counter for operation number (1 to 48)
    #sym_count = 1 # count number of symmetries

    #? ============================================== Rotations ========================================================
    #print("============== Rotations ===============")
    for axis, theta in zip(rot_axes, rot_angles):
        k = axis / np.linalg.norm(axis) # axis of rotation
        rot_vertices = list([rodrigues_rotation(np.array((x, y, z)), k, theta) for (x, y, z) in second_shape])
        shifted_coord = shift_coordinates(rot_vertices)
        #print(shifted_coord)
        if sorted_first_shape == sorted(shifted_coord):
            #sym_count+=1
            #print(f"i={operation_num[i]}", '\u2714', axis, np.around(math.degrees(theta),1))
            return 1 # the shapes are same 
        #else: print(f"i={operation_num[i]}", '\u2717', axis, np.around(math.degrees(theta),1))
        #i+=1

    #? ============================================= Reflections =======================================================
    #print("============== Reflections =============")
    for normal, d in zip(ref_planes, d_list):
        normal = normal / np.linalg.norm(normal) # normal of plane
        ref_vertices = list([reflect_across_plane(np.array((x, y, z)), normal, d) for (x, y, z) in second_shape])
        shifted_coord = shift_coordinates(ref_vertices)
        shifted_coord = list([tuple(np.around((x,y,z), 4)) for (x,y,z) in shifted_coord])
        #print(shifted_coord)

        if sorted_first_shape == sorted(shifted_coord):
            #sym_count+=1
            #print(f"i={operation_num[i]}", '\u2714', np.around(normal,4), np.around(d,4))
            return 1  
        #else: print(f"i={operation_num[i]}", '\u2717', np.around(normal,4), np.around(d,4))
        #i+=1

    #? ============================================ Rotoreflections ====================================================
    #print("============== H = D x K ===========")
    # Operation H = D x K
    vertex_axes = rot_axes[15:]
    vertex_angles = rot_angles[15:]
    for axis, theta in zip(vertex_axes, vertex_angles):
        k = axis / np.linalg.norm(axis) # axis of rotation
        rot_vertices = list([rodrigues_rotation(np.array((x, y, z)), k, theta) for (x, y, z) in second_shape])
        inverted_coord = inversion(rot_vertices)
        shifted_coord = shift_coordinates(inverted_coord)
        #print(shifted_coord)

        if sorted_first_shape == sorted(shifted_coord): 
            #sym_count+=1
            #print(f"i={operation_num[i]}", '\u2714', axis, np.around(math.degrees(theta),1),"D x K")
            return 1
        #else: print(f"i={operation_num[i]}", '\u2717', axis, np.around(math.degrees(theta),1), "D x K")
        #i+=1

    #? Operation J = A x E =========================================================================================
    #print("============== J = A x E ===========")
    for j in range(len(rotoJ_axes)):
        # rotation A
        axis = rotoJ_axes[j]; theta = rotoJ_angles[j]
        k = axis / np.linalg.norm(axis) # axis of rotation
        rot_vertices = list([rodrigues_rotation(np.array((x, y, z)), k, theta) for (x, y, z) in second_shape])
        
        # roto reflection A x E
        normal = rotoJ_planes[j]; d = rotoJ_d_list[j]
        normal = normal / np.linalg.norm(normal) # normal of plane #! Not really necessary here
        ref_vertices = list([reflect_across_plane(np.array((x, y, z)), normal, d) for (x, y, z) in rot_vertices])
        shifted_coord = shift_coordinates(ref_vertices)
        shifted_coord = list([tuple(np.around((x,y,z), 4)) for (x,y,z) in shifted_coord])
        #print(shifted_coord)
        if sorted_first_shape == sorted(shifted_coord): 
            #sym_count+=1
            #print(f"i={operation_num[i]}", '\u2714',  np.around(math.degrees(theta),1), np.around(normal,4), np.around(d,4), "A x E")
            return 1
        #else: print(f"i={operation_num[i]}", '\u2717',  np.around(math.degrees(theta),1), np.around(normal,4), np.#around(d,4), "A x E")
        #i+=1

    #? Inversion ====================================================================================================
    ref_vertices = inversion(second_shape)
    shifted_coord = shift_coordinates(ref_vertices)
    shifted_coord = list([tuple(np.around((x,y,z), 4)) for (x,y,z) in shifted_coord])
    #print(shifted_coord)
    if sorted_first_shape == sorted(shifted_coord): 
        #sym_count+=1
        #print(f"i={operation_num[i]}", '\u2714', np.around(normal,4), np.around(d,4), "K")
        return 1
    #else: print(f"i={operation_num[i]}", '\u2717', np.around(normal,4), np.around(d,4), "K")
    #print("Order =", sym_count)
    return 0

#========================================================================================================#
def compare_polycubes(tile_coord, valid_shape_i):
    """ 
    Checks if the newly assembled shape is similar to an already encountered shape.
    This executes some elimination checks to speed up comparison and then does brute_force_comparison.

    Args:
        - tile_coord (list of tuples) : tile coordinates of the new polycube.
        - valid_shape_i (list of tuples) : ith 'valid' shape to compare with the new shape.

    Returns:
        - 1 if the new shape is equivalent to the ith valid shape under any symmetry operation.
        - 0 otherwise.

    Raises:
        - None
    """
    A = list(valid_shape_i) # list of tile coordinates of ith valid shape
    B = list(tile_coord) # list of tile coordinates of new encountered shape
    #print("A \n{} \nB \n{}".format(A,B))

    if len(A) != len(B): return 0 # unequal number of tiles = dissimilar shapes
    
    else: # number of tiles in A = B
        #! INSPECT lengths_A and lengths_B are correctly calculated
        A_new = shift_coordinates(A) # all non-negative coordinates of ith valid shape
        lengths_A = return_lengths(A)[0] # lengths of A
        B_new = shift_coordinates(B) # all non-negative coordinates of newly encountered shape
        lengths_B = return_lengths(B)[0] # lengths of B

        if set(A_new) == set(B_new): return 1 # identical shape, exit
        if sorted(lengths_A) != sorted(lengths_B): return 0 # range along 3 axes don't match, exit
            
        output = brute_force_comparison(tile_coord, valid_shape_i) 
        return output
    