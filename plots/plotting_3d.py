"""
# Description: Contains functions for plotting cubes of a polycube.

# Functions:
    - plot_all_cubes(params, ax, coord_list, picked_tiles, cube_outline, axes_lines)
    - plot_all_polycubes(params, path, filepath, sort)

# Dependencies:
    - numpy
    - matplotlib.pyplot
    - matplotlib.colors
    - mpl_toolkits.mplot3d.art3d
    - utils.return_lengths
"""

# Import inbuilt packages
import numpy as np
import matplotlib 
#matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection

# Import user defined packages
from utils import return_lengths, convert_tilecoord_to_2d_matrix

#====================================================================================================================#
def plot_shape_mat(ax3,shape_matrix):
    """
    For 2D polyominoes, this plots the shape matrix.

    Args:
        - ax3 (axes): axes to plot the shape matrix
        - shape_matrix (2D array): shape matrix of the polyomino

    Returns:
        - 2D plot of the shape matrix
    """
    #fig, ax3 = plt.subplots()
    n_rows = len(shape_matrix); n_cols = len(shape_matrix[0])

    # Plot matrix
    cmap = plt.get_cmap('Blues')#.copy()
    cmap.set_under('white')
    #shape_mat_flip = np.flipud(shape_matrix) #? NEW ADDITION: it is necessary only when shape_mat has been flipped once
    ax3.matshow(shape_matrix, cmap=cmap, vmin = 0.5, origin='lower')

    # Minor ticks
    ax3.set_xticks(np.arange(-.5, n_cols, 1), minor=True)
    ax3.set_yticks(np.arange(-.5, n_rows, 1), minor=True)

    # Gridlines based on minor ticks
    ax3.grid(which='minor', color='w', linestyle='-', linewidth=2)

    # Remove minor ticks
    ax3.tick_params(which='minor', bottom=False, left=False, top=False)
    ax3.tick_params(which='major', bottom=False, left=False, top=False)
    ax3.set_xticklabels([]); ax3.set_yticklabels([])

    plt.box(False)

#====================================================================================================================#
def plot_all_cubes(params, ax, coord_list, picked_tiles=[], cube_outline='False', axes_lines='False', dim=3):
    """ 
    Plots all cubes.

    Version: 1.0 [29 Jan 2026]
    Author: Prarthana Agrawal

    Older versions:
        v0.1
        - Added dim parameter to allow plotting of 2D shape matrices and 3D cubes using the same function.

    Args:
        - params (dict): {'tile_types':tile_types, 'side_types':side_types, 'neutral_sides':neutral_sides, 'self_int_sides':self_int_sides, 'n_rules':n_rules, 'dim':dim}
            tile_types : number of tile species
            side_types : number of side types
            neutral_sides : number of neutral sides
            self_int_sides : number of self-interacting sides
            n_rules : number of rules
            dim : dimension of the polycube assembly
        - ax (axes3D)
        - coord_list (list of tuples) : origins of all assembled cubes
        - picked_tiles (list): list of names of all placed cubes in order of being assembled
        - cube_outline (str): to display cube outline 'True' or 'False'
        - axes_lines (str): to display X,Y,Z axes lines 'True' or 'False'
        - dim (int): dimension of the shape (2 or 3)

    Returns:
        - 2D plot of shape matrix (for dim=2), or
        - 3D plot of all cubes with a lightsource

    Example use:
        for 3D plotting:
            %matplotlib widget
            fig = plt.figure()
            ax1 = fig.add_subplot(111, projection='3d')
            ax1.view_init(elev=-149, azim=138)
        
        for 2D plotting:
            fig2, ax2 = plt.subplots()
    """

    #============================================= Initializations =================================================#
    if len(picked_tiles) == 0:
        picked_tiles = ['A']*len(coord_list)

    # Unpack parameters
    if params == {}: # if no params are passed, assume only 1 tile type
        tile_types = ['A']
    else:
        tile_types = params['tile_types']
    n_tiles = len(tile_types)

    #===============================================================================================================#
    # ---------------------------------------- For 2D polyominoes plotting ------------------------------------------
    #===============================================================================================================#
    if dim == 2: 
        #fig, ax = plt.subplots()
        shape_matrix = convert_tilecoord_to_2d_matrix(coord_list)
        plot_shape_mat(ax,shape_matrix)
        return
    
    #===============================================================================================================#
    # ---------------------------------------- For 3D polycubes plotting -------------------------------------------
    #===============================================================================================================#

    # --------------------------------------------- Color of tiles -------------------------------------------------- #
    cmap = plt.get_cmap('Blues').copy()
    colors = cmap(np.linspace(0.1,0.8,n_tiles))
    tile_color_dict = dict()
    for tl in range(n_tiles):
        tile_color_dict[tile_types[tl]] = colors[tl]
    light = LightSource(azdeg=270, altdeg=30) # Create a Light Source
    
    lengths, min_vals, max_vals = return_lengths(coord_list)

    #----------------------------------------- Loop over each tile/cube ----------------------------------------------#
    for tile_num, (x0,y0,z0) in enumerate(coord_list):
       
        # Coordinates of vertices
        shr = 0.2 # shrink to create some gaps between cubes (better visualisation)
        x = [x0+a for a in [0,1-shr,0,0,1-shr,0,1-shr,1-shr]]
        y = [y0+a for a in [0,0,1-shr,0,1-shr,1-shr,0,1-shr]]
        z = [z0+a for a in [0,0,0,1-shr,0,1-shr,1-shr,1-shr]]

        # Triangular connections
        vertices = [[0,1,4],[0,4,2],[2,3,5],[2,0,3],[1,6,4],[6,7,4],[0,1,3],[1,3,6],[3,5,7],[3,6,7]]  
        tupleList = list(zip(x, y, z))
        poly3d = [[tupleList[vertices[ix][iy]] for iy in range(len(vertices[0]))] for ix in range(len(vertices))]
        ax.add_collection3d(Poly3DCollection(np.asarray(poly3d), edgecolors=tile_color_dict[picked_tiles[tile_num]], facecolors=tile_color_dict[picked_tiles[tile_num]], linewidths=1, alpha=1, zorder=2, shade=True, lightsource=light))
        ax.set_box_aspect([1,1,1])  

        # Outline of each cube
        if cube_outline == 'True':
            ax.add_collection3d(Line3DCollection([
            [(x0+a, y0+b, z0+c) for a, b, c in [(0,0,0), (1-shr,0,0), (1-shr,1-shr,0), (0,1-shr,0), (0,0,0)]],  # Bottom square
            [(x0+a, y0+b, z0+c) for a, b, c in [(0,0,1-shr), (1-shr,0,1-shr), (1-shr,1-shr,1-shr), (0,1-shr,1-shr), (0,0,1-shr)]],  # Top square
            [(x0+a, y0+b, z0+c) for a, b, c in [(0,0,0), (0,0,1-shr)]],  # Vertical lines
            [(x0+a, y0+b, z0+c) for a, b, c in [(1-shr,0,0), (1-shr,0,1-shr)]],
            [(x0+a, y0+b, z0+c) for a, b, c in [(1-shr,1-shr,0), (1-shr,1-shr,1-shr)]],
            [(x0+a, y0+b, z0+c) for a, b, c in [(0,1-shr,0), (0,1-shr,1-shr)]]],
            colors='lightgray', linewidths=1,zorder=3))
            
    #---------------------------------------- Drawing X, Y, Z axes lines ---------------------------------------------#
    if axes_lines == 'True':
        ext = 2 # measure of length of axes
        (max_x, max_y, max_z) = max_vals
        (min_x, min_y, min_z) = min_vals
        pos_axes = [[(0,0,0),(max_x+ext,0,0)],[(0,0,0),(0,max_y+ext,0)],[(0,0,0),(0,0,max_z+ext)]]
        neg_axes = [[(0,0,0),(min_x-ext,0,0)],[(0,0,0),(0,min_y-ext,0)],[(0,0,0),(0,0,min_z-ext)]]
        ax.add_collection3d(Line3DCollection(pos_axes, colors='gray', linewidths=1, zorder=1))
        ax.add_collection3d(Line3DCollection(neg_axes, colors='gray', linewidths=1, linestyles='dashed', zorder=1))
        ax.text(max_x+ext,0,0,"X", fontsize=20, color='gray')
        ax.text(0,max_y+ext,0,"Y", fontsize=20, color='gray')
        ax.text(0,0,max_z+ext,"Z", fontsize=20, color='gray')
    
    #---------------------------------------------- Setting axes limits --------------------------------------------#
    minval = min(min_vals); maxval = max(max_vals)
    ax.set_xlim(minval, maxval+1); ax.set_ylim(minval, maxval+1); ax.set_zlim(minval, maxval+1)

    plt.axis('off'); plt.grid(b=None)
    plt.tight_layout();
    #plt.show()

    
#====================================================================================================================#
def plot_all_polycubes(params, path, filepath, sort, dim=3, first_n_shapes=None, target_size=None, target_symorder=None, target_symclass=None, print_shapes=False):
    """
    Plot all polycubes in a grid format after reading the data from the shape files.
    This can plot them in decreasing order of frequency.

    Version: 1.0 [06 Feb 2026]

    Older versions:
        v0.3
        - Added dim parameter to allow plotting of 2D shape matrices and 3D cubes using the same function.
        - Added target_symorder and target_symclass parameters to filter shapes based on symmetry.
        - Added print_shapes parameter to optionally print shape coordinates to the console.

        v0.2 [24 Feb 2025]
        - changed xlimit calculations.
        - added option to plot first n number of shapes. 

    Args:
        - params (dict): {'tile_types':tile_types, 'side_types':side_types, 'neutral_sides':neutral_sides,  'self_int_sides':self_int_sides, 'n_rules':n_rules, 'dim':dim}
            tile_types : number of tile species
            side_types : number of side types
            neutral_sides : number of neutral sides
            self_int_sides : number of self-interacting sides
            n_rules : number of rules
            dim : dimension of the polycube assembly
        - path (str): path to the data files
        - filepath (str): file path to the data files
        - sort (bool): sort the shapes based on frequency True or False
        - dim (int): dimension of the shape (2 or 3)
        - first_n_shapes (int): number of shapes to plot (eg first 64 shapes)
        - target_size (int): size of the polycube to plot (say focus on only 16-mers)
        - target_symorder (int): target symmetry order to filter shapes
        - target_symclass (str): target symmetry class to filter shapes
        - print_shapes (bool): whether to print shape coordinates to the console

    Returns:
        - 2D plot of all polyominoes' shape matrices in a grid format, or
        - 3D plot of all polycubes in a grid format    

    Raises:
        - ValueError: If frequency and shape data do not match
        - ValueError: If no shapes are found
    """
    # ============================================= Import packages ================================================== #
    import math
    from more_itertools import sort_together 
    from matplotlib.gridspec import GridSpec 
    from utils import shift_coordinates 
    from symmetry import get_lunnon_data, get_classes_and_orders_dict 
    from analyzer.compactness import calculate_radius_of_gyration 

    # ============================================= Load data ======================================================== #
    frequency = np.loadtxt(path + filepath + 'frequency.txt')
    complexity = np.loadtxt(path + filepath + 'complexity.txt')
    complexity_species = np.loadtxt(path + filepath + 'complexity_species.txt')
    lz_complexity = np.loadtxt(path + filepath + 'lz_complexity.txt')
    shape_file = open(path + filepath + 'valid_shapes.txt', 'r')
    shape_content = shape_file.readlines()

    try: 
        symmetry_classes = np.loadtxt(path + filepath + 'symmetry_classes.txt', dtype=str)
        symmetry_orders = np.loadtxt(path + filepath + 'symmetry_orders.txt', dtype=int)
    except FileNotFoundError:
        print("Symmetry class and order files not found. Assigning default values of 'X' for symmetry class and 0 for symmetry order.")
        symmetry_classes = ['X'] * len(shape_content)  # Create a list of 'X' with the same length as shape_content
        symmetry_orders = [0] * len(shape_content) # No symmetry classes available

    # ============================================ Data validation ================================================== #
    # Check if frequency and shape data match
    if len(frequency) != len(shape_content): raise ValueError("Error: Frequency and shape data do not match.")
    # Check if shapes are found
    if len(frequency) == 0: raise ValueError("Error: No shapes found.")

    # Read shapes, shift coordinates then store in a list
    shapes = []
    for line in shape_content:
        coords = eval(line.strip())
        shift_coords = shift_coordinates(coords)
        shapes.append(shift_coords)

    # Sort shapes based on frequency
    if sort:
        freq_sort, shape_sort, symmetry_classes_sort, symmetry_orders_sort, complexity_sort, complexity_species_sort, lz_complexity_sort = sort_together([frequency, shapes, symmetry_classes, symmetry_orders, complexity, complexity_species, lz_complexity], reverse=True)
    else:
        freq_sort = frequency
        shape_sort = shapes
        symmetry_classes_sort = symmetry_classes
        symmetry_orders_sort = symmetry_orders
        complexity_sort = complexity
        complexity_species_sort = complexity_species
        lz_complexity_sort = lz_complexity

    # ============================================ Filter shapes =================================================== #
    def filter_criteria(shape, shape_index, size_criteria=False, sym_order_criteria=False, sym_class_criteria=False):
        """
        Filter criteria based on which parameters are provided (True/False).
        E.g. if size = True, filter based on size.
             if sym_order = True, filter based on symmetry order.
             if sym_class = True, filter based on symmetry class.
        
        Args:
            - shape (list of tuples): coordinates of the polycube
            - size_criteria (bool): filter based on size
            - sym_order_criteria (bool): filter based on symmetry order
            - sym_class_criteria (bool): filter based on symmetry class

        Returns:
            - bool: True if the shape meets the filter criteria, False otherwise
        """
        meets_criteria = True

        if size_criteria:
            meets_criteria &= (len(shape) == target_size)
        if sym_order_criteria:
            meets_criteria &= (symmetry_orders_sort[shape_index] == target_symorder)
        if sym_class_criteria:
            meets_criteria &= (symmetry_classes_sort[shape_index] == target_symclass)
        return meets_criteria

    size_criteria = False; sym_order_criteria = False; sym_class_criteria = False
    if target_size is not None: size_criteria = True
    if target_symorder is not None: sym_order_criteria = True
    if target_symclass is not None: sym_class_criteria = True
    
    # Filter shapes based on the target size
    if size_criteria or sym_order_criteria or sym_class_criteria:

        filtered_shapes = []
        filtered_freqs = []
        filtered_symmetry_classes = []
        filtered_symmetry_orders = []
        filtered_complexity = []
        filtered_complexity_species = []
        filtered_lz_complexity = []
        shape_index = 0

        for shape, freq, sym_class, sym_order, comp, comp_species, lz_comp in zip(shape_sort, freq_sort, symmetry_classes_sort, symmetry_orders_sort, complexity_sort, complexity_species_sort, lz_complexity_sort):

            if filter_criteria(shape, shape_index, size_criteria, sym_order_criteria, sym_class_criteria):
                filtered_shapes.append(shape)
                filtered_freqs.append(freq)
                filtered_symmetry_classes.append(sym_class)
                filtered_symmetry_orders.append(sym_order)
                filtered_complexity.append(comp)
                filtered_complexity_species.append(comp_species)
                filtered_lz_complexity.append(lz_comp)

            shape_index += 1

        if len(filtered_shapes) == 0:
            raise ValueError("Error: No shapes found matching the filter criteria.")
        
        shape_sort = filtered_shapes
        freq_sort = filtered_freqs
        symmetry_classes_sort = filtered_symmetry_classes
        symmetry_orders_sort = filtered_symmetry_orders
        complexity_sort = filtered_complexity
        complexity_species_sort = filtered_complexity_species
        lz_complexity_sort = filtered_lz_complexity

    # ------------------------------------------------------------------------------------------------------------- #
    # ============================================ Plotting ======================================================= #
    # ------------------------------------------------------------------------------------------------------------- #
    # Number of shapes to plot
    minval = []; maxval = []
    if first_n_shapes is not None:
        n_shapes = min(first_n_shapes, len(shape_sort))
    else: 
        n_shapes = len(shape_sort)

    # Get min and max values for axes limits
    max_vals = dict()
    for shape in shape_sort[:n_shapes+1]:
        lengths = return_lengths(shape)[0]
        len_x, len_y, len_z = lengths
        max_vals['x'] = max(len_x, max_vals.get('x', 0))
        max_vals['y'] = max(len_y, max_vals.get('y', 0))
        max_vals['z'] = max(len_z, max_vals.get('z', 0))
    maxval = max(max_vals['x'], max_vals['y'], max_vals['z'])

    # Determine number of rows and columns for subplots
    if n_shapes > 6 and n_shapes < 8:
        cols = 4
        rows = 2
    else:
        cols = math.ceil(math.sqrt(n_shapes))
        rows = math.ceil(n_shapes / cols)
    
    fig = plt.figure(figsize=(cols * 3, rows * 3))  # Adjust the figure size as needed
    gs = GridSpec(rows, cols, figure=fig)

    # Create a supertitle
    if params == {}:
        pass
    else:
        n_tiles = len(params['tile_types'])
        n_sides = len(params['side_types'])
        dim = params['dim']
        n_rules = params['n_rules']
        if n_rules == 'all':
            n_rules_str = n_rules
        else:
            n_rules_str = f"$10^{int(np.log10(n_rules))}$"
        supertitle = f"{dim}D {n_tiles}s{n_sides}c {n_rules_str}"
        fig.suptitle(supertitle, fontsize=20)

    # Plotting all shapes
    for i in range(n_shapes):

        tile_coords = shape_sort[i] # coordinates of the polycube
        if print_shapes: print(tile_coords)
        
        row = i // cols # row number
        col = i % cols # column number
        
        # calculate radius of gyration
        # Rg = calculate_radius_of_gyration(tile_coords)

        if dim == 3:
            ax = fig.add_subplot(gs[row, col], projection='3d') # add subplot
            ax.view_init(elev=-149, azim=138) # set view angle
        elif dim == 2:
            ax = fig.add_subplot(gs[row, col]) # add subplot

        if dim == 2:
            ax.set_xlim(0, maxval)
            ax.set_ylim(-1, maxval)
            ax.set_box_aspect(1)
            ax.set_aspect('equal')

        plot_all_cubes(params, ax, tile_coords, picked_tiles=[], cube_outline='False', axes_lines='False', dim=dim)

        #ax.set_xlim(0, max_vals['x']); ax.set_ylim(0, max_vals['y']); ax.set_zlim(0, max_vals['z'])
        if dim == 3: 
            ax.set_xlim(0, maxval); ax.set_ylim(0, maxval); ax.set_zlim(0, maxval)
        
        #--------------------------- Set title with complexity and frequency ------------------------------#
        freq_percent = np.around((freq_sort[i]/np.sum(freq_sort))*100, 2) # frequency percentage

        if symmetry_classes_sort[i] == 'X':
            title = f"{int(complexity_species_sort[i])}$\\tilde{{K}}_s$ {int(complexity_sort[i])}$\\tilde{{K}}_c$ {int(lz_complexity_sort[i])}$\\tilde{{K}}_{{lz}}$ \nf={freq_percent}%"
            
        else:
            symmetry_class_label = symmetry_classes_sort[i]
            symmetry_order_label = symmetry_orders_sort[i]
            
            if dim == 2:
                title = (
                    f"{int(complexity_species_sort[i])}$\\tilde{{K}}_s$ "
                    f"{int(complexity_sort[i])}$\\tilde{{K}}_c$ "
                    f"{int(lz_complexity_sort[i])}$\\tilde{{K}}_{{lz}}$\n"
                    f"f={freq_percent}%, "
                    f"$\\mathrm{{{symmetry_class_label}}}$"
                )
            else:
                title = (
                    f"{int(complexity_species_sort[i])}$\\tilde{{K}}_s$ "
                    f"{int(complexity_sort[i])}$\\tilde{{K}}_c$ "
                    f"{int(lz_complexity_sort[i])}$\\tilde{{K}}_{{lz}}$\n"
                    f"f={freq_percent}%, "
                    f"$\\mathrm{{{symmetry_class_label}}}_{{{symmetry_order_label}}}$"
                )
            
        if dim==3: ax.set_title(title, fontsize=14, loc='right')
        if dim==2: ax.set_title(title, fontsize=12, loc='center', y=0.9)

        #ax.axis('on')

    # Adjust spacing between subplots to control gaps (tweak wspace/hspace as needed)
    #fig.subplots_adjust(wspace=0.05, hspace=0.05,bottom=0.05)
    plt.show()
    #return ax