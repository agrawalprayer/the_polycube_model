"""
Utilities for binary operations.

[Author: Prarthana Agrawal]
[Date: 11 Feb 2025]
[Version: 0.1]

Functions:
    - rqrd_num_bits(n_sides): Determine the number of bits needed to represent a set of numbers in binary.
    - number_to_binary(n, num_bits): Convert a number to binary representation.
    - orientation_to_binary(orientation): Convert orientation to binary representation.
    - tiledict_orientdict_to_binary(n_sides, tile_dict, orient_dict): Convert a genotype (tile_dict, orient_dict) to binary representation.
    - genotype_lists_to_binary(nbits, flat_genotype, flat_orientation): Convert a genotype list and orientations to a concatenated binary genotype string.
    - lempel_ziv_complexity78(sequence): Calculate the Lempel-Ziv complexity for a binary sequence. 1978
    - lempel_ziv_complexity76(sequence): Calculate the Lempel-Ziv complexity using the 1976 parsing algorithm.
    - lempel_ziv_complexity76_symmetric(x): Calculate the symmetric Lempel-Ziv complexity for a binary sequence using the 1976 parsing algorithm.
"""

# Import packages
import math

# Functions
def rqrd_num_bits(n_sides):
    """ 
    Determine the number of bits needed to represent a set of numbers in binary.
    The numbers go from 0,1,2,..., n_sides-1

    Args:
        - n_sides (int): The maximum value in the set of numbers (starting from 0)

    Returns:
        - int: The number of bits needed to represent the numbers in binary

    """
    return math.ceil(math.log2(n_sides))  # Optimal bits needed

#================================================================================================#
def number_to_binary(n, num_bits):
    """ 
    Convert a number to binary representation.

    Args:
        - n (int): The number to convert
        - num_bits (int): The number of bits to use for the binary representation
    
    Returns:
        - str: The binary representation of the number
    """
    return format(n, f'0{num_bits}b')

#================================================================================================#
def orientation_to_binary(orientation):
    """ 
    Convert orientation to binary representation.

    Args:
        - orientation (str): The orientation to convert (valid inputs are 'O', 'U', 'I', 'D', 'L', 'R')

    Returns:
        - str: The binary representation of the orientation

    Raises:
        - ValueError: If the orientation is not valid
    """
    if orientation not in ['O', 'U', 'I', 'D', 'L', 'R']:
        raise ValueError(f"Invalid orientation: {orientation}")
    
    orientation_to_binary_dict = {
    'O': '000',
    'U': '001',
    'I': '010',
    'D': '011',
    'L': '100',
    'R': '101'
    }

    #Note: An alternate way to achieve this is to use 1 bit for sign and 2 bits to represent directions  LR, UD, IO
    return orientation_to_binary_dict[orientation]

#================================================================================================#
def tiledict_orientdict_to_binary(n_sides, tile_dict, orient_dict, dim):
    """ 
    Convert a genotype (tile_dict, orient_dict) to binary representation.

    Args:
        - n_sides (int): The maximum value in the set of numbers (starting from 0)
        - tile_dict (dict): The genotype dictionary.
        - orient_dict (dict): The orientation dictionary.
    
    Returns:
        - str: The binary representation of the genotype.
    """
    # Determine the number of bits needed to represent each patch
    num_bits = rqrd_num_bits(n_sides)
    #print(f'Number of bits per patch: {num_bits}')
    binary_genotype = '' # Initialize the binary genotype (combined tiledict and orientdict)

    if dim not in (2, 3):
        raise ValueError("dim must be 2 or 3")
    
    # Iterate over each tile
    for key in tile_dict.keys():
        interface = tile_dict[key]
        if dim == 3:
            orientation = orient_dict[key]

        # Iterate over each interface and orientation
        for i in range(len(interface)):  # should be 6
            binary_interface = number_to_binary(int(interface[i]), num_bits)
            if dim==3: binary_orientation = orientation_to_binary(orientation[i])
            # Combine binary interface and orientation
            if dim == 3: binary_patch = str(binary_interface) + str(binary_orientation)
            else: binary_patch = str(binary_interface)
            #print(f'{key}: {interface[i]} -> {binary_interface} -> {binary_orientation} -> {binary_patch}')
            # Append the binary patch to the binary genotype
            binary_genotype = binary_genotype + binary_patch
            
    return binary_genotype

#================================================================================================#
def genotype_lists_to_binary(nbits, flat_genotype, flat_orientation, dim):
    """ 
    Convert a genotype list (flattened, containing all Nt tiles) and orientations (flattened, containing all Nt
    tiles) to a concatenated binary genotype string.

    Args:
        - nbits (int): Number of bits required to represent integers to binary.
        - flat_genotype (list): flattened genotype list
        - flat_orientation (list): flattened orientation list
        - dim (int): The dimension of the genotype (e.g., 3 for 3D)
    
    Returns:
        - str: The binary representation of the genotype.
    """
    if dim not in (2, 3):
        raise ValueError("dim must be 2 or 3")
    
    if dim == 3:
        if flat_orientation is None:
            raise ValueError("flat_orientation must be provided for 3D")
        if len(flat_orientation) != len(flat_genotype):
            raise ValueError("Genotype and orientation lengths differ")

    binary_genotype = '' # Initialize the binary genotype (combined tiledict and orientdict)

    # One by one convert every element of genotype and orientations to binary and concat
    for i in range(len(flat_genotype)):
        
        # Load interface and orientation values
        interface_patch = flat_genotype[i]
        if dim==3: orient_patch = flat_orientation[i] if flat_orientation is not None else None
        #print(f'Interface: {interface_patch}, Orientation: {orient_patch}')

        # Convert interface and orientation values to binary
        binary_interface_patch = number_to_binary(int(interface_patch), nbits)
        if dim==3: binary_orient_patch = orientation_to_binary(orient_patch) if orient_patch is not None else ''
        #print(f'Binary interface: {binary_interface_patch}, Binary orientation: {binary_orient_patch}')

        # Combine binary interface and orientation
        if dim==3: binary_patch = str(binary_interface_patch) + str(binary_orient_patch)
        else: binary_patch = str(binary_interface_patch)
        #print(f'Binary patch: {binary_patch}')

        # Append the binary patch to the binary genotype
        binary_genotype = binary_genotype + binary_patch
        
    return binary_genotype

#================================================================================================#
def lempel_ziv_complexity78(sequence):
    """
    #! Note this is closer to LZ78 implementation than LZ76.
    #! Deprecated 13 Jul 2026.
    Calculate the Lempel-Ziv complexity for a binary sequence.

    Args:
        sequence (str): The binary sequence to compute complexity for. It also works for non-binary sequences.

    Returns:
        int: The Lempel-Ziv complexity of the input sequence.

    Description:
        The function iterates through the sequence, extracting substrings starting at the current position with increasing length.
        If the substring is already in the set of seen substrings, the length is incremented to consider a longer substring.
        If the substring is not in the set, it is added to the set, and the starting position is updated to move past the current substring.

    Adapted from: 
        https://github.com/Naereen/Lempel-Ziv_Complexity/blob/master/Short_study_of_the_Lempel-Ziv_complexity.ipynb

    Example usage:
        sequency = '00100100110111000'
        lz_complexity = lempel_ziv_complexity(sequency)
        print(f"LZ complexity: {lz_complexity}")
    """
    
    seen_substrings = set()
    sequence_length = len(sequence)
    start_index = 0
    substring_length = 1
    complexity_count = 0

    while True:
        # Break the loop if the end of the sequence is reached
        if start_index + substring_length > sequence_length:
            break

        # Extract the current substring
        current_substring = sequence[start_index : start_index + substring_length]

        # If the substring is already seen, increase the length to consider a longer substring
        if current_substring in seen_substrings:
            substring_length += 1
        else:
            # Add the new substring to the set of seen substrings
            seen_substrings.add(current_substring)
            # Move the start index past the current substring
            start_index += substring_length
            # Reset the substring length to 1
            substring_length = 1
            # Increment the complexity count
            complexity_count += 1

    return complexity_count

#================================================================================================#
def lempel_ziv_complexity76(s):
    """
    N_w(x): number of words (distinct patterns) in the dictionary produced by the Lempel-Ziv (1976) parsing algorithm, computed via the Kaspar & Schuster (1987) two-pointer method. Works on any indexable sequence (str, tuple, list), not just binary strings.
    Date: 13 Jul 2026
    """
    n = len(s)
    if n == 0:
        return 0
    if n == 1:
        return 1
    i, k, l = 0, 1, 1
    c = 1
    k_max = 1
    stop = False
    while not stop:
        if s[i + k - 1] != s[l + k - 1]:
            if k > k_max:
                k_max = k
            i += 1
            if i == l:
                c += 1
                l += k_max
                if l + 1 > n:
                    stop = True
                else:
                    i, k, k_max = 0, 1, 1
            else:
                k = 1
        else:
            k += 1
            if l + k > n:
                c += 1
                stop = True
    return c


def lempel_ziv_complexity76_symmetric(x):
    """
    C_LZ(x): the Lempel-Ziv-based complexity estimator.
    Strictly binary: x must be a string over {'0','1'}.

        C_LZ(x) = log2(n)                                    if x = 0^n or 1^n
                = log2(n) * (N_w(x) + N_w(reverse(x))) / 2    otherwise

    where N_w is the number of words (distinct patterns) produced by the Lempel-Ziv (1976) parsing (Kaspar & Schuster 1987 two-pointer algorithm), and n = |x|.
    Date: 13 Jul 2026
    """
    import math

    n = len(x)
    if n == 0:
        return 0.0  # empty sequence (e.g. unassembled single-cube state): complexity 0
    if not set(x) <= {"0", "1"}:
        raise ValueError("x must be a binary string over {'0','1'}")
    if n == 1:
        return 0.0  # log2(1) = 0, matches the x=0^n / 1^n branch for n=1

    if x == "0" * n or x == "1" * n:
        return math.log2(n)
    fwd = lempel_ziv_complexity76(x)
    rev = lempel_ziv_complexity76(x[::-1])
    return math.log2(n) * (fwd + rev) / 2.0