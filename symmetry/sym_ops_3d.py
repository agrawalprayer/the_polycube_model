"""
# Description: functions for symmetry operations in 3D 
    - The symmetry operations are 24 rotations and 24 reflections+rotoreflections of a cube.

# Functions:
    - rodrigues_rotation(v, k, theta): Rotate vector v by angle theta(radian) around axis k using Rodrigues' rotation formula.
    - reflect_across_plane(point, normal, d): Reflects a point across a plane defined by its normal vector and distance from the origin.
    - inversion(coord_list): Inverts points under bbox's center (x,y,z) --> (-x,-y,-z)

# Dependencies:
    - numpy
"""

# Import packages
import numpy as np

def rodrigues_rotation(v, k, theta):
    """
    Rotate vector v by angle theta(radian) around axis k using Rodrigues' rotation formula.
    
    Args:
        - v (np.array): The vector to be rotated.
        - k (np.array): The unit vector representing the axis of rotation.
        - theta (float): The angle of rotation in radians.
        
    Returns:
        - np.array: The rotated vector.
    """
    v = np.array(v)
    k = np.array(k)
    k = k / np.linalg.norm(k)  # Ensure k is a unit vector

    v_rot = (v * np.cos(theta) +
             np.cross(k, v) * np.sin(theta) +
             k * np.dot(k, v) * (1 - np.cos(theta)))
    
    return np.round(v_rot,3) #! Note the rounding off here! #! is this necessary?

#==============================================================================================================#
def reflect_across_plane(point, normal, d):
    """
    Reflects a point across a plane defined by its normal vector and distance from the origin.
    
    Args:
        - point (array): The point to reflect (numpy array).
        - normal (array): The normal vector of the plane (numpy array).
        - d (float): The plane constant (distance from origin).
        
    Returns: 
        - The reflected point (tuple).
    """
    normal = normal / np.linalg.norm(normal)  # Normalize the normal vector
    distance = (np.dot(normal, point) - d)
    reflected_point = point - 2 * distance * normal
    return tuple(reflected_point)

#==============================================================================================================#
def inversion(coord_list):
    """ 
    Inverts points under bounding box's center (x,y,z) --> (-x,-y,-z)

    Args:
        - coord_list (list of tuples): cube origins of the polycube

    Returns:
        - list of tuples: inverted cube origins
    """
    inverted_coord = [(-x,-y,-z) for (x,y,z) in coord_list]
    return inverted_coord


