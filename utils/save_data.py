"""
# Description: Contains functions for saving data files generated in polycube self-assembly.

# Functions:
    - save_all_genotypes(file_name, all_genotypes, all_orientations, save_orientations=False)
    - save_assembly_description(file_name, params, assembly_settings, parallel_run_config, master_seed, sol_stats)

# Dependencies:
    - None
"""

def save_all_genotypes(file_name, all_genotypes, all_orientations, save_orientations=False):
    """ 
    To save all genotypes and orientations that have been assembled.

    Args:
        - file_name (str): name of the data file.
        - all_genotypes (list): list of all genotypes.
        - all_orientations (list): list of all orientations.
        - save_orientations (bool): whether to save orientations or not.

    Returns:
        - None
    """
    # Store all genotypes
    with open(file_name + 'all_genotypes.txt', 'w') as f:
        for genotype in all_genotypes:
            for element in genotype:
                f.write("{} ".format(element))
            f.write("\n")

    # Store all orientations -- meaningless for 2D
    if save_orientations:
        with open(file_name + 'all_orientations.txt', 'w') as f:
            for orient in all_orientations:
                for element in orient:
                    f.write("{} ".format(element))
                f.write("\n")

#======================================================================================================================#
def save_assembly_description(file_name, params, assembly_settings, parallel_run_config, master_seed, sol_stats):
    """
    To save details about the assembly in a data file.

    Args:
        - file_name (str): name of the data file
        - params (dict): {'tile_types':tile_types, 'side_types':side_types, 'neutral_sides':neutral_sides, 'self_int_sides':self_int_sides, 'n_rules':n_rules, 'dim':dim}
        - assembly_settings (dict): {'assembly_type': assembly_type, 'Dmax': Dmax, 'max_tiles': max_tiles, 'kmax': kmax}
        - parallel_run_config (dict): {'tot_splits': tot_splits, 'split': split}
        - sol_stats (dict): {'UBD': UBD, 'ND': ND, 'valid': valid}
            UBD : number of unbounded solutions
            ND : number of non-deterministic solutions
            valid : number of valid solutions
        - master_seed (int): seed for random number generator

    Returns:
        - None
    """
    # load parameter values
    tile_types = params['tile_types']
    side_types = params['side_types']
    n_tiles = len(tile_types)
    n_sides = len(side_types)

    neutral_sides = params['neutral_sides']
    self_int_sides = params['self_int_sides']
    n_rules = params['n_rules']
    dim = params['dim']

    assembly_type = assembly_settings['assembly_type']
    Dmax = assembly_settings['Dmax']
    kmax = assembly_settings['kmax']
    max_tiles = assembly_settings['max_tiles']

    tot_splits = parallel_run_config['tot_splits']
    #nsplit = parallel_run_config['nsplit']
    
    UBD = sol_stats['UBD']
    ND = sol_stats['ND']
    valid = sol_stats['valid']

    file = open(file_name + 'description.txt', 'w')
    file.writelines(['*'*50])
    file.writelines(['\nNumber of tile types = ', str(n_tiles)])
    file.writelines(['\nNumber of side types = ', str(n_sides)])
    file.writelines(['\nNeutral sides = ', str(neutral_sides)])
    file.writelines(['\nSelf interacting sides = ', str(self_int_sides)])
    file.writelines(['\n'+'*'*50])
    file.writelines(['\nDimensions = ', str(dim)])
    file.writelines(['\nNumber of rules checked = ', str(n_rules)])
    file.writelines(['\n'+'*'*50])
    file.writelines(['\nAssembly type = ', str(assembly_type)])
    file.writelines(['\nDmax = ', str(Dmax)])
    file.writelines(['\nmax_tiles = ', str(max_tiles)])
    file.writelines(['\nkmax = ', str(kmax)])
    file.writelines(['\nmaster_seed = ', str(master_seed)])
    file.writelines(['\ntot_splits = ', str(tot_splits)])
    file.writelines(['\n'+'*'*50])
    file.writelines(['\nNumber of unbounded rules = ', str(UBD)])
    file.writelines(['\nNumber of non deterministic rules = ', str(ND)])
    file.writelines(['\nNumber of valid rules = ', str(valid)])
    file.writelines(['\n'+'*'*50])
    
    file.writelines(['\n'])
    file.close()