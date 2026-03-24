"""
# Description: Contains key assembly-setup functions like tile identities and attachment rules.

# Functions:
    - import_input(input_file_path): Import the input.py file from the specified path.
    - get_params(n_tiles, n_sides, neutral_sides, self_int_sides, n_rules, dim): Returns system parameters.
    - tile_identity_func(params): Generate tile identities with randomly assigned interfaces and orientations.
    - rules_func(params): Generate rules for side interactions based on given criteria.

# Dependencies:
    - random
    - importlib.util
    - os

"""

# Import packages
import random

def import_input(input_file_path):
    """
    Function to import the input.py file from the specified path.

    Args:
        - input_file_path (str): Path to the input.py file.

    Returns:
        - input_module: Imported module from the input.py file.

    Raises:
        - FileNotFoundError: If the input file does not exist.    
    """

    import os # required for file operations
    import importlib.util # required for dynamic loading of input.py

    # Check if the file exists
    if not os.path.isfile(input_file_path):
        print(f"Error: {input_file_path} does not exist.")
        return None

    # Dynamically load the input.py file from the specified path
    spec = importlib.util.spec_from_file_location("input", input_file_path)
    input_module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(input_module)
    
    # Return the imported module
    return input_module

#=============================================================================================#
def get_params(n_tiles, n_sides, neutral_sides, self_int_sides, n_rules, dim):
    """
    Returns system parameters for self-assembly process in a two-digit format '00'.

    Args:
        - n_tiles (int) : number of tile types / species
        - n_sides (int) : number of side types / colors
        - neutral_sides (list) : neutral sides (single-digit-format eg [0,3])
        - self_int_sides (list) : self-interacting sides (single-digit-format eg [3])
        - n_rules (int) : number of genotypes to assemble
        - dim (int=2/3) : dimensionality of assembly, 2D or 3D    
    
    Returns:
        - params (dict) : {'tile_types':tile_types, 'side_types':side_types, 'neutral_sides':neutral_sides,         'self_int_sides':self_int_sides, 'n_rules':n_rules, 'dim':dim}
    """
    
    tile_types = [chr(x) for x in range(65,65+n_tiles)] # types of tiles
    side_types = list(range(0,n_sides)) # types of sides/interfaces "Always take even number of sides excluding 0"
    formatter = "{:02d}".format
    side_types = list(map(formatter,side_types))
    neutral_sides = list(map(formatter, neutral_sides))
    self_int_sides = list(map(formatter, self_int_sides))

    #print("Tile types{} \nSide types {} \nNeutral sides {} \nSelf-interacting sides  {}".format(tile_types, side_types, neutral_sides, self_int_sides)) 
    #print("="*50)

    params = {
        'tile_types': tile_types,
        'side_types': side_types,
        'neutral_sides': neutral_sides,
        'self_int_sides': self_int_sides,
        'n_rules': n_rules,
        'dim': dim
    }
    return params

#=============================================================================================#
def tile_identity_func(params):
    """
    Generate tile identities with randomly assigned interfaces and orientations.
    
    Args:
        - params (dict) : {'tile_types':tile_types, 'side_types':side_types, 'neutral_sides':neutral_sides,         'self_int_sides':self_int_sides, 'n_rules':n_rules, 'dim':dim}
            tile_types (list): List of tile names eg ['A', 'B'].
            side_types (list): List of interface types eg ['00', '01', '02'].
            neutral_sides (list): List of neutral sides eg ['00', '03'].
            self_int_sides (list): List of self-interacting sides eg ['03'].
            n_rules (int): Number of rules to generate.
            dim (int): Dimension of assembly (2 or 3).
        
    Returns:
        - tile_dict (dict): Dictionary mapping tile names to their interface types.
        - orient_dict (dict): Dictionary mapping tile names to their orientations.

    Raises:
        - ValueError: If invalid value of dim is provided.
    """
    tile_types = params['tile_types']
    side_types = params['side_types']
    dim = params['dim']

    tile_dict = dict() # key:tile name, value:interfaces in North-East-South-West-Back-Front format 
    orient_dict = dict() # key:tile name, value:string of six letters indicating directions of arrows

    for tile in tile_types:

        # 3D Assembly ---------------------------------------------------------------------------------
        if dim == 3:

            # Generate random side types for 3D tiles
            sides = random.choices(side_types, k=6) # randomly assign six interface-types to each cube
            tile_dict[(tile)] = sides

            # Generate random orientations for 3D tiles In(I), Out(O), Left(L), Right(R), Up(U), Down(D)
            North = random.choices(['I','O','L','R'], k=1)[0]
            East  = random.choices(['I','O','U','D'], k=1)[0]
            South = random.choices(['I','O','L','R'], k=1)[0]
            West  = random.choices(['I','O','U','D'], k=1)[0]
            Back  = random.choices(['U','D','L','R'], k=1)[0]
            Front = random.choices(['U','D','L','R'], k=1)[0]
            orient_dict[(tile)] = North + East + South + West + Back + Front

        # 2D Assembly ---------------------------------------------------------------------------------
        elif dim == 2:
            sides = random.choices(side_types, k=4) # randomly assign NESW interfaces
            tile_dict[(tile)] = sides + ['00','00'] # neutralise the back and front patches for 2D assembly
            orient_dict[(tile)] = 'OOOOUU' # All NESW outwards, last two don't matter as they are neutralised
            
        else: raise ValueError("Invalid value of dim (Enter 2 or 3)")

    #print("Tiles and their interfaces are", tile_dict)
    #print("Tiles and their orientations are", orient_dict)
    return tile_dict, orient_dict

#=============================================================================================#
def rules_func(params):
    """ 
    Generate rules for side interactions based on given criteria.
        - odd integer i connects with even integer i+1 eg 1 <--> 2, 3 <--> 4.
        - neutral sides do not connect with anything.
        - self-interacting side connects with itself.

    Args:
        - params (dict) : {'tile_types':tile_types, 'side_types':side_types, 'neutral_sides':neutral_sides,         'self_int_sides':self_int_sides, 'n_rules':n_rules, 'dim':dim}
            tile_types (list): list of all tiles.
            side_types (list): list of all sides.
            neutral_sides (list): list of sides that do not connect.
            self_int_sides (list): list of sides that connect with themselves.
            n_rules (int): number of rules to generate.
            dim (int): dimension of assembly (2 or 3).

    Returns:
        - rules_dict (dict): dictionary mapping each side to its interaction partner eg {1:2, 2:1, 3:4, 4:3, 5:5}.

    Raises:
        - Exception: If number of non-self interacting sides is not even.
    """

    side_types = params['side_types']
    neutral_sides = params['neutral_sides']
    self_int_sides = params['self_int_sides']
    rules_dict = dict()

    # Filter out neutral and self-interacting sides
    nonself_int_sides = list(
        filter(lambda side: (side not in neutral_sides and side not in self_int_sides), side_types))
    #print("Non-self interacting sides ", nonself_int_sides)
  
    # Raise error if it is impossible to form unique pairs of interacting sides
    if len(nonself_int_sides)%2 != 0:
        raise Exception("Error: number of non-self interacting sides is not even to form unique pairs")
    
    # Store connections eg 1 <--> 2, 3 <--> 4, 5 <--> 6 
    for indx in range(0, len(nonself_int_sides), 2):
        rules_dict[nonself_int_sides[indx]] = nonself_int_sides[indx+1]
        rules_dict[nonself_int_sides[indx+1]] = nonself_int_sides[indx]

    # Store any self-interactions if present eg 7 <--> 7
    if len(self_int_sides) > 0:
        for side in self_int_sides:
            rules_dict[side] = side

    #print("Rules dict ", rules_dict)
    return rules_dict
