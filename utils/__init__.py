"""
# Description: This __init__.py file is used to treat the 'utils' directory as a Python package.
    - It imports specific functions or classes from the modules within the package to expose them at the package level.
    - This allows for easier and cleaner imports when using the package.

# Modules imported:
    - polycube_utils.py
    - save_data.py
    - print_utils.py

# Public API:
The __all__ list defines the public API of the package, specifying which functions or classes are available for import when using 'from utils import *'.
 
These can then be imported as
    from utils import return_length
instead of
    from utils.polycube_utils import return_length
"""

from .polycube_utils import *
from .save_data import *
from .print_utils import *
#from .binary_utils import *
from .read_data import *

# Define the public API using __all__
__all__ = [
    "return_length",
    "extract_underscore",
    "shift_coordinates",
    "get_bounding_box",
    "get_batch_range",
    "get_exh_genotypes",
    "combine_parallel_runs",
    "save_all_genotypes",
    "save_assembly_description",
    "print_msg_box",
    "extract_sol_stats",
    "get_bar_shift",
    "zero_singleton_frequencies",
    "convert_tilecoord_to_2d_matrix",]
    #"rqrd_num_bits",
    #"number_to_binary",
    #"orientation_to_binary",
    #"genotype_to_binary",
    #"lempel_ziv_complexity",
#]