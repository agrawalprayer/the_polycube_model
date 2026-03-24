
"""
# Description: This __init__.py file is used to treat the 'symmetry' directory as a Python package.
    - It imports specific functions or classes from the modules within the package to expose them at the package level.
    - This allows for easier and cleaner imports when using the package.

# Modules imported:
    - sym_ops_3d.py
    - compare_shapes.py

# Public API:
The __all__ list defines the public API of the package, specifying which functions or classes are available for import when using 'from symmetry import *'.
"""

from .sym_ops_3d import *
from .compare_shapes import *
from .symmetry_analysis import *

# Define the public API using __all__
__all__ = [
    "rodrigues_rotation",
    "reflect_across_plane",
    "inversion",
    "brute_force_comparison",
    "compare_polycubes",
    "permitted_operations",
    "get_symmetry_order",
    "get_2d_symmetry_order_and_class",
    "save_symmetry_orders",
    "get_symmetry_class",
    "get_lunnon_data",
    "sym_classes_per_order",
    "get_classes_and_orders_dict",
    "get_symclass_to_vector_dict",
    "convert_class_to_vector"
]