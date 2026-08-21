"""
# Description: Contains all integer-sequences from OEIS database related to polycubes (2D and 3D).

# Functions:
    - get_oeis_data(dim, type_of_polycube)
    - get_oeis_symmetry_data(dim, type_of_polycube, sym_class)
    - get_oeis_3Dclass_matrix()
    - get_oeis_3Dorder_matrix()

# Dependencies:
    - numpy

"""

# Import packages
import numpy as np


def get_oeis_data(dim, type_of_polycube):
    """ 
    Return oeis data about number of polycubes of a certain size for given dimension and type.

    Args:
        - dim (int): Dimension of polycube
        - type_of_polycube (str): Type of polycube (e.g. 'free', 'fixed', 'onesided')
    
    Returns:
        - polycube_size_n (list): List of polycube sizes
        - num_of_polycubes_of_size_n (list): List of number of polycubes of size n
        - seq_name (str): Name of the sequence in OEIS database

    """
    #==================================================================================================#
    #----------------------------------------  2D Polyominoes -----------------------------------------#
    #==================================================================================================#
    #! Certain sequences start with n=0 (which doesn't make sense physically so n=0 is dropped here).

    if dim == 2:

        if type_of_polycube == 'free':

            # https://oeis.org/A000105/list
            seq_name = 'A000105'
            num_of_polycubes_of_size_n = [1,1,2,5,12,35,108,369,1285,4655,17073,63600,
            238591,901971,3426576,13079255,50107909,192622052,
            742624232,2870671950,11123060678,43191857688,
            168047007728,654999700403,2557227044764,
            9999088822075,39153010938487,153511100594603]

        elif type_of_polycube == 'onesided':

            # https://oeis.org/A000988/list
            seq_name = 'A000988'
            num_of_polycubes_of_size_n = [1,1,2,7,18,60,196,704,2500,9189,33896,126759,
            476270,1802312,6849777,26152418,100203194,
            385221143,1485200848,5741256764,22245940545,
            86383382827,336093325058,1309998125640,
            5114451441106,19998172734786,78306011677182,
            307022182222506,1205243866707468,4736694001644862]

        elif type_of_polycube == 'fixed':

            #https://oeis.org/A001168/list
            seq_name = 'A001168'
            num_of_polycubes_of_size_n = [1,2,6,19,63,216,760,2725,9910,36446,135268,
            505861,1903890,7204874,27394666,104592937,
            400795844,1540820542,5940738676,22964779660,
            88983512783,345532572678,1344372335524,
            5239988770268,20457802016011,79992676367108,
            313224032098244,1228088671826973]
            
        else: raise ValueError("Invalid type of polycube. Use 'free', 'onesided' or 'fixed'.")
            
    #==================================================================================================#
    #----------------------------------------  3D Polycubes -----------------------------------------#
    #==================================================================================================#
    elif dim == 3:

        if type_of_polycube == 'free':

            # https://oeis.org/A038119/list
            seq_name = 'A038119'
            num_of_polycubes_of_size_n = [1,1,2,7,23,112,607,3811,25413,178083,1279537,
            9371094,69513546,520878101,3934285874,29915913663,
            228779330204,1758309223457,13573319825615,
            105192814197984,818136047201932,6383528588447574]

        elif type_of_polycube == 'onesided':

            # https://oeis.org/A000162/list
            seq_name = 'A000162'
            num_of_polycubes_of_size_n = [1,1,2,8,29,166,1023,6922,48311,346543,2522522,
            18598427,138462649,1039496297,7859514470,
            59795121480,457409613979,3516009200564,
            27144143923583,210375361379518,1636229771639924,
            12766882202755783]

        elif type_of_polycube == 'fixed':

            # https://oeis.org/A001931/list
            seq_name = 'A001931'
            num_of_polycubes_of_size_n = [1,3,15,86,534,3481,23502,162913,1152870,8294738,
            60494549,446205905,3322769321,24946773111,
            188625900446,1435074454755,10977812452428,
            84384157287999,651459315795897,5049008190434659,
            39269513463794006,306405169166373418]

        else: raise ValueError("Invalid type of polycube. Use 'free', 'onesided' or 'fixed'.")

    else: raise ValueError("Invalid dimension. Use 2 or 3.")

    polycube_size_n = np.arange(1,len(num_of_polycubes_of_size_n)+1,1)

    return polycube_size_n, num_of_polycubes_of_size_n, seq_name


#==================================================================================================#

def get_oeis_3d_symmetry_data(sym_class):
    """ 
    Return oeis data about number of polycubes having size n and symmetry class from Lunnon's 33 classes of symmetry.
    
    Args:
        - sym_class (str): from Lunnon's 33 classes of symmetry

    Returns:
        - n_list (list): List of polycube sizes
        - num_of_polycubes_of_size_n (list): List of number of polycubes of size n
        - seq_name (str): Name of the sequence in OEIS database

    # free polycubes https://oeis.org/A038119 (allowing reflections)
    """
    
    if sym_class == 'I':
        seq_name = 'A376964' # https://oeis.org/A376964
        num_of_polycubes_of_size_n = [0,0,0,0,4,46,394,3025,22707,167732,1241417,
        9221624,68936674,518574100,3925132946,29878869619,
        228629549175,1757697391087,13570818452472,
        105182527335313,818093680980786,6383353461322488]

    elif sym_class == 'A':
        seq_name = 'A376965' # https://oeis.org/A376965
        num_of_polycubes_of_size_n = [0,0,0,0,0,0,0,0,0,0,0,2,0,2,2,22,14,31,41,213,182,321,453,1796]

    elif sym_class == 'B':
        seq_name = 'A376966'
        num_of_polycubes_of_size_n = [0,0,0,0,0,3,4,37,52,342,502,2836,4343,22622,35405,176176,281141,1363112,2205171,10527712, 17221126,81462884,134424679,632308448]

    elif sym_class == 'C':
        seq_name = 'A376967'
        num_of_polycubes_of_size_n = [0,0,0,1,2,5,17,49,138,374,1062,2851,8010,21432, 60142,161386,453253,1222110,3436564,9316409, 26231463,71496106]

    elif sym_class == 'D':
        seq_name = 'A376968'
        num_of_polycubes_of_size_n = [0,0,0,0,0,0,1,0,1,9,1,7,68,9,61,473,81,440,3316,614,3315,23537,4613,24419,169404,34247,182074, 1236200,255636,1363691,9131384]

    elif sym_class == 'E':
        seq_name = 'A376969'
        num_of_polycubes_of_size_n = [0,0,0,1,6,28,126,520,2099,8429,33676,135201,543248,2195182,8893547,36196788,147739046, 605138811,2485051070,10233284681,42234214910,174693725151]

    elif sym_class == 'F':
        seq_name = 'A376970'
        num_of_polycubes_of_size_n = [0,0,0,0,3,6,26,72,241,623,2028,5374,16781,45094,138971,378196,1156028,3183390,9679303,26929415, 81575661,228958790]

    elif sym_class == 'H':
        seq_name = 'A376972'
        num_of_polycubes_of_size_n = [0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,10,4,6,1,0,0,78,31,47,11,3,0,565,212]

    elif sym_class == 'J':
        seq_name = 'A376973'
        num_of_polycubes_of_size_n = [0,0,0,0,0,0,0,0,0,1,1,6,3,10,12,63,41,97,124,553,401,821,1097,4464]
    
    elif sym_class == 'K':
        seq_name = 'A376974'
        num_of_polycubes_of_size_n = [0,0,0,0,0,1,1,20,17,192,175,1632,1522,13088,12339,102846,97374,803363,760594,6274638,5928747, 49133524]
    
    elif sym_class == 'BB':
        seq_name = 'A376975'
        num_of_polycubes_of_size_n = [0,0,0,0,0,0,0,0,0,1,1,5,3,13,16,58,47,136,181,546,521,1247,1729,4698]
    
    elif sym_class == 'BC':
        seq_name = 'A376976'
        num_of_polycubes_of_size_n = [0,0,0,0,0,0,0,0,0,1,2,7,5,17,23,78,62,181,218,717,591,1606,1835,6002]
    
    elif sym_class == 'BE':
        seq_name = 'A376977'
        num_of_polycubes_of_size_n = [0,0,0,1,1,5,5,21,26,93,116,392,497,1616,2087,6690,8656,27493,35797,113728,147838,469824,611551,1951237]

    elif sym_class == 'BF':
        seq_name = 'A376978'
        num_of_polycubes_of_size_n = [0,0,0,0,0,1,0,2,2,12,9,37,32,134,114,426,375,1393,1236,4398,3934,13883]

    elif sym_class == 'CE':
        seq_name = 'A376979'
        num_of_polycubes_of_size_n = [0,0,1,0,3,4,12,13,48,61,186,231,727,941,2855,3707,11293,14909,44955,59760,180092,241827]

    elif sym_class == 'CK':
        seq_name = 'A376980'
        num_of_polycubes_of_size_n = [0,0,0,0,0,2,1,6,5,26,17,69,50,235,159,693,470,2142,1399,6350,4086,19273]

    elif sym_class == 'EE':
        seq_name = 'A376981'
        num_of_polycubes_of_size_n = [0,0,0,1,2,6,11,32,61,158,308,749,1481,3481,6997, 16086,32692,74114,152114,341410,705978,1573062, 3273527,7253791]

    elif sym_class == 'CD':
        seq_name = 'A376982'
        num_of_polycubes_of_size_n = [0,0,0,0,0,0,0,0,0,1,0,1,0,1,1,4,2,10,5,7,8,20,13,81,39,60,65,123,98,598,269]

    elif sym_class == 'FF':
        seq_name = 'A376983'
        num_of_polycubes_of_size_n = [0,0,0,1,0,0,2,0,1,6,2,3,15,6,14,46,22,48,133,66, 161,395,218,522,1162,673,1651,3465,2105,5279,10369]

    elif sym_class == 'AB':
        seq_name = 'A377127'
        num_of_polycubes_of_size_n = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,3,0,0,0,14,4,2]

    elif sym_class == 'AE':
        seq_name = 'A376984'
        num_of_polycubes_of_size_n = [0,0,0,0,0,0,0,1,0,0,0,3,2,0,2,15,9,0,9,57,40,2,40,238]

    elif sym_class == 'BFF':
        seq_name = 'A376985'
        num_of_polycubes_of_size_n = [0,0,0,0,0,0,0,1,0,0,0,4,1,1,2,14,6,5,9,52,22,21]

    elif sym_class == 'CJ':
        seq_name = 'A376986'
        num_of_polycubes_of_size_n = [0,0,0,0,0,1,1,2,3,5,6,11,16,25,33,56,81,118,162,263,380,556,768,1229]

    elif sym_class == 'EEE':
        seq_name = 'A376987'
        num_of_polycubes_of_size_n = [0,0,0,0,0,1,2,3,4,7,12,16,29,39,73,92,170,209,405,485,932,1095,2166,2528]
    
    elif sym_class == 'EF':
        seq_name = 'A376988'
        num_of_polycubes_of_size_n = [0,0,0,0,0,1,1,2,3,6,7,13,20,32,43,71,109,162,226,363,553,811]
    
    elif sym_class == 'EFF':
        seq_name = 'A376989'
        num_of_polycubes_of_size_n = [0,0,0,0,0,0,1,2,1,2,4,9,8,13,19,42,37,65,78,181,152,298]

    elif sym_class == 'BD':
        seq_name = 'A377128'
        num_of_polycubes_of_size_n = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,1,0,1,0,0,4,2,2,3,4,0,23,7,5,10,3,1,48,16,19,28,49,2,174,49,58,84,46,18,406,111,169,238,424,34,1285,321,524,678,410,153,3139,747,1393,1872,3185]

    elif sym_class == 'CF':
        seq_name = 'A376990'
        num_of_polycubes_of_size_n = [0,0,0,0,0,1,0,0,0,0,0,3,2,2,1,0,0,8,4,6,4,1,0,24,11,15,10,4,1,77,32]

    elif sym_class == 'BBC':
        seq_name = 'A376991'
        num_of_polycubes_of_size_n = [0,1,1,2,2,1,1,2,4,2,5,7,8,5,10,15,20,12,23,36,48,30]

    elif sym_class == 'CCC':
        seq_name = 'A377129'
        num_of_polycubes_of_size_n = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,1,0,0,0,0,2,1,0,0,3,0,2,2,2,0,1,1,9,2,0,0,14,1,7,5,10,1,4,4,31,6,6,4,42,4,25,13,45,9,15,13,111,20,28,21,143,14,95,44,175,34,64,44,401,68,111,76,482]

    elif sym_class == 'DEE':
        seq_name = 'A377130'
        num_of_polycubes_of_size_n = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,3,1,0,0,0,1,6,2,1,0,0,4,13,5,1,0,0,6,28,9,4,0,0,20,61,26,7,0,0,36,129,43,18,0,0,94,274,109,33,0,0,182,582,201,81,2,0,438,1231,501]

    elif sym_class == 'R':
        seq_name = 'A377131'
        num_of_polycubes_of_size_n = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,3,0,0,0,0,0,0,0,12]

    elif sym_class == 'G':
        seq_name = 'A376971'
        num_of_polycubes_of_size_n = [1,0,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,2,1,0,0,0,1,2,1,1,0,0,1,2,3,1,0,0,1,3,1,1,0,0,4,5,4,1,0,0,6,7,4,3,0,0,8,10,11,3,0,0,12,14,8,5,1,0,22,21,21,7,0,0,34,32,20,12,2,0,50,48,48,16,1,1,76,69,48,27,8,1]

    else: raise ValueError("Invalid symmetry class. Use one of Lunnon's 33 classes of symmetry.")

    n_list = np.arange(1,len(num_of_polycubes_of_size_n)+1,1) 

    return n_list, num_of_polycubes_of_size_n, seq_name


#==================================================================================================#
def get_oeis_3Dclass_matrix():
    """ 
    Return a matrix of Lunnon's symmetry classes with number of polycubes of size n for each class.
    Rows = symmetry classes, Columns = polycube sizes
    Entries = number of polycubes of size n for each symmetry class.

    Returns:
        - oeis_class_matrix (np.array): Matrix of number of polycubes of size n for each symmetry class    
    """
    # List of all Lunnon's 33 classes of symmetry with their orders
    from symmetry import get_lunnon_data
    lunnon_classes_orders_dict = get_lunnon_data(sort=True)
    lunnon_sym_classes_list = list(lunnon_classes_orders_dict.keys())

    # First find what is the smallest value of n for which all symmetry classes data is available
    min_n = 100
    for i, sym_class in enumerate(lunnon_sym_classes_list):
        n_list, num_of_polycubes_of_size_n, seq_name = get_oeis_3d_symmetry_data(sym_class)
        if max(n_list) < min_n:
            min_n = max(n_list)
    #print("minimum n for which all symmetry classes have available oeis data is", min_n)

    # Create a matrix with rows as symmetry classes and columns as number of polycubes of size n (min avaliable n) 
    oeis_class_matrix = np.zeros((len(lunnon_sym_classes_list), min_n), dtype=np.int64) #(33x22)

    for i, sym_class in enumerate(lunnon_sym_classes_list):
        n_list, num_of_polycubes_of_size_n, seq_name = get_oeis_3d_symmetry_data(sym_class)
        oeis_class_matrix[i][:min_n] = num_of_polycubes_of_size_n[:min_n]
    #print(oeis_matrix)

    return oeis_class_matrix

#==================================================================================================#
def get_oeis_3Dorder_matrix():
    """ 
    Return a matrix of 3D symmetry orders with number of polycubes of size n for each order.
    Rows = symmetry orders, Columns = polycube sizes
    Entries = number of polycubes of size n for each symmetry order.

    Returns:
        - oeis_order_matrix (np.array): Matrix of number of polycubes of size n for each symmetry order    
    
    Note: OEIS data is available for each symmetry 'class' and not 'order'. So first, we combine classes with same order. Then we combine the data for each class to get the data for each order.
    """
    from collections import defaultdict

    # List of all Lunnon's 33 classes of symmetry with their orders
    from symmetry import get_lunnon_data
    lunnon_classes_orders_dict = get_lunnon_data(sort=True)
    lunnon_sym_classes_list = list(lunnon_classes_orders_dict.keys())

    # -----------------------------------------------------------------------------------------------#
    # ------------------------- Combine the classes with same order ---------------------------------#
    # -----------------------------------------------------------------------------------------------#
    # Create a defaultdict of lists
    grouped_order_sym_class_dict = defaultdict(list)

    # Iterate through the original dictionary
    for key, value in lunnon_classes_orders_dict.items():
        grouped_order_sym_class_dict[value].append(key)

    # Convert defaultdict to a regular dictionary and sort the keys
    grouped_order_sym_class_dict = dict(sorted(grouped_order_sym_class_dict.items()))


    # -----------------------------------------------------------------------------------------------#
    # ------------------------- Combine the data for each order -------------------------------------#
    # -----------------------------------------------------------------------------------------------#
    oeis_order_matrix = np.zeros((10,100), dtype=np.int64) # order x n (upto size 100)

    for i, (sym_order, sym_classes) in enumerate(grouped_order_sym_class_dict.items()):

        combined_n_list = None
        combined_num_of_polycubes_list = None

        for sym_class in sym_classes:
            n_list, num_of_polycubes_of_size_n, seq_name = get_oeis_3d_symmetry_data(sym_class)
            
            if combined_n_list is None:
                combined_n_list = n_list  # Start with the first available list
                combined_num_of_polycubes_list = np.array(num_of_polycubes_of_size_n)
            else:
                # for each order, identify the smallest available n_list
                min_length = min(len(combined_n_list), len(n_list))
                
                # Truncate to smallest available n_list so far
                combined_n_list = combined_n_list[:min_length]
                combined_num_of_polycubes_list = combined_num_of_polycubes_list[:min_length] + np.array(num_of_polycubes_of_size_n[:min_length])

        combined_num_of_polycubes_list = list(combined_num_of_polycubes_list)
        oeis_order_matrix[i][:max(combined_n_list)] = combined_num_of_polycubes_list

    return oeis_order_matrix


#==================================================================================================#
def get_oeis_2D_symmetry_data(sym_class):
    """ 
    Return oeis data about number of polyominoes having size n and symmetry class from {C1, C2, C4, D1, D2, D4}
    
    Args:
        - sym_class (str): from 6 classes of symmetry for 2D polyominoes  

    Returns:
        - n_list (list): List of polyomino sizes
        - num_of_polyominoes_of_size_n (list): List of number of polyominoes of size n
        - seq_name (str): Name of the sequence in OEIS database

    """
    
    if sym_class == 'D4':
        seq_name = "A142886"
        num_of_polyominoes_of_size_n = [1,0,0,1,1,0,0,1,2,0,0,3,2,0,0,5,4,0,0,12,7,0,0,20,11,0,0,45,20,0,0,80,36,0,0,173,65,0,0,310,117,0,0,664,216,0,0,1210,396,0,0,2570,736,0,0,4728, 1369,0,0,9976,2558,0,0,18468,4787,0,0,38840]

    elif sym_class == 'D2':
        seq_A = "A056877"
        part_A = [0,1,1,1,1,2,3,4,4,8,10,15,17,30,35,60,64,117,128,
        236,241,459,476,937,912,1813,1789,3706,3456,7187,
        6779,14712,13161,28571,25839,58457,50348,113798,
        98957,232718,193375,453969,380522,927601,745248,
        1813219,1468202,3702063]

        seq_B = "A056878"
        part_B = [0,0,0,0,0,0,1,1,0,1,2,3,3,5,6,14,9,20,20,56,32,
        80,64,224,114,315,217,863,397,1234,751,3331,1400,
        4816,2632,12815,4973,18792,9349,49400,17810,73338,
        33557,190643,64309,286368,121511,737532,233891,
        1119215,443271,2859154]

        seq_name = seq_A + ',' + seq_B
        num_of_polyominoes_of_size_n = [a + b for a, b in zip(part_A, part_B)]

    elif sym_class == 'D1':
        seq_A = "A006746"
        part_A = [0,0,0,1,2,6,9,23,38,90,147,341,564,1294,2148,
        4896,8195,18612,31349,70983,120357,271921,463712,
        1045559,1792582,4034832,6950579,15619507,27023509,
        60638559,105320716,236006955,411364068,920626423,
        1609836928]

        seq_B = "A006748"
        part_B = [0,0,1,0,2,2,7,5,26,22,91,79,326,301,1186,1117,
        4352,4212,16119,15849,60174,60089,226146,228426,
        854803,872404,3247207,3342579,12389106,12850662,
        47448984,49544820,182338754,191529007,702807040,
        742163178,2716205709,2882119756]
        
        seq_name = seq_A + ',' + seq_B
        num_of_polyominoes_of_size_n = [a + b for a, b in zip(part_A, part_B)]

    elif sym_class == 'C4':
        seq_name = "A144553"
        num_of_polyominoes_of_size_n = [0,0,0,0,0,0,0,1,0,0,0,3,2,0,0,12,7,0,0,44,25,0,0,165,90,0,0,603,319,0,0,2235,1136,0,0,8283,4088,0, 0,30936,14868,0,0,116111,54526,0,0,438465,201527,0,0,1663720,750169,0,0,6342211,2809931,0,0, 24273767]

    elif sym_class == 'C2':
        seq_name = "A006747"
        num_of_polyominoes_of_size_n = [0,0,0,1,1,5,4,18,19,73,73,278,283,1076,1090,4125,
        4183,15939,16105,61628,62170,239388,240907,932230,
        936447,3641945,3651618,14262540,14277519,55987858,
        55961118,220223982,219813564,867835023,865091976,
        3425442681]

    elif sym_class == 'C1':
        # subtract total - all other symmetries to get the number of polyominoes with no symmetry (C1)
        #seq_name = "A006479" #! this is wrongly referenced in symmetry and simplicity paper
        n_tot, a_n_tot, seq_tot = get_oeis_data(dim=2, type_of_polycube='free')
        n_D4, a_n_D4, seq_D4 = get_oeis_2D_symmetry_data(sym_class='D4')
        n_D2, a_n_D2, seq_D2 = get_oeis_2D_symmetry_data(sym_class='D2')
        n_D1, a_n_D1, seq_D1 = get_oeis_2D_symmetry_data(sym_class='D1')
        n_C4, a_n_C4, seq_C4 = get_oeis_2D_symmetry_data(sym_class='C4')
        n_C2, a_n_C2, seq_C2 = get_oeis_2D_symmetry_data(sym_class='C2')

        seq_name = "A006479 = " + seq_tot + " - " + seq_D4 + " - " + seq_D2 + " - " + seq_D1 + " - " + seq_C4 + " - " + seq_C2
        num_of_polyominoes_of_size_n = [tot - D4 - D2 - D1 - C4 - C2 for tot, D4, D2, D1, C4, C2 in zip(a_n_tot, a_n_D4, a_n_D2, a_n_D1, a_n_C4, a_n_C2)]
        
    else: raise ValueError("Invalid symmetry class. Use one of {C1, C2, C4, D1, D2, D4}.")

    n_list = list(np.arange(1,len(num_of_polyominoes_of_size_n)+1,1))
    return n_list, num_of_polyominoes_of_size_n, seq_name