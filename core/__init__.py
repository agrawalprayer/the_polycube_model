"""
# Description: This __init__.py file is used to treat the 'core' directory as a Python package.
    - It imports specific functions or classes from the modules within the package to expose them at the package level.
    - This allows for easier and cleaner imports when using the package.

# Modules imported:
    - system_setup.py
    - assembly_helpers.py

# Public API:
The __all__ list defines the public API of the package, specifying which functions or classes are available for import when using 'from core import *'.
"""

from .system_setup import *
from .assembly_helpers import *

# Define the public API using __all__
__all__ = [
    "import_input",
    "get_params",
    "tile_identity_func",
    "rules_func",
    "origin_finder",
    "zero_sides",
    "path_finder"]
