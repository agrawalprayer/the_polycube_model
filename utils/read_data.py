"""
# Description: Contains functions for reading and extracting values from saved data files.

# Functions:
    - extract_sol_stats(path, n_tiles, n_sides): Extract solution stat data from description files.

# Dependencies:
    - re # for regular expressions
"""

def extract_sol_stats(path, n_tiles, n_sides):
    """ 
    Open the description files and return solution stat data.

    Args:
        - path (str): path to the directory containing the main folder.
        - n_tiles (int): number of tile types.
        - n_sides (int): number of sides per tile.

    Returns:
        - sol_stats (dict): dictionary containing the number of UBD, ND and valid solutions.

    [Date: 24 Feb 2025]
    [Author:Prarthana Agrawal]
    """
    import re
    
    filename = f'data_files/{n_tiles}s{n_sides}c_combined_description.txt'
    with open(path+filename, 'r') as file:
        data = file.read()
    
    # Regular expressions to extract values
    ubd_match = re.search(r'(?:Total number|Number) of unbounded rules = (\d+)', data)
    nd_match = re.search(r'(?:Total number|Number) of non deterministic rules = (\d+)', data)
    valid_match = re.search(r'(?:Total number|Number) of valid rules = (\d+)', data)

    # Extract values if matches are found
    ubd = int(ubd_match.group(1)) if ubd_match else None
    nd = int(nd_match.group(1)) if nd_match else None
    valid = int(valid_match.group(1)) if valid_match else None
    
    # Store values in a dictionary
    sol_stats = {'UBD': ubd, 'ND': nd, 'Valid': valid}
    return sol_stats