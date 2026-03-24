"""
# Description: This __init__.py file is used to treat the 'plots' directory as a Python package.
    - It imports specific functions or classes from the modules within the package to expose them at the package level.
    - This allows for easier and cleaner imports when using the package.

# Modules imported:
    - plotting_3d.py

# Public API:
The __all__ list defines the public API of the package, specifying which functions or classes are available for import when using 'from plots import *'.
"""

from .plotting_3d import *

# Define the public API using __all__
__all__ = [
    "plot_all_cubes",
    "plot_all_polycubes"
]