# the_polycube_model
A computational model of polycube self-assembly, in which cubes connect via specified interaction rules to form discrete structures. The model also applies to 2D polyomino assembly.

## How it works:
1. Start by specifying the number of different tile types (`n_tiles`). These are labelled as 'A', 'B', 'C'... in the input dictionary. 
    Note: This means current implementation would throw an error beyond n_tiles > 26. But this can be easily fixed by changing how input dictionaries are handled.
2. Specify the number of different interface types (`n_sides`). These are stored as string of two digit numbers '00', '01' ,'02'.. in the input dictionary. 
    Note: This implies that this is limited to maximum number '99'. This can also be easily extended further by not restricting to this format.
3. Specify which interface types are 'neutral' or 'self-interacting'.
4. Fix the dimension of self-assembly: 2D or 3D.


## How 3D assembler works for 2D:
Each tile has six sides by default. To switch to 2D mode of assembly, the last two sides are considered '00' (neutral). This refers to 'back' and 'front' faces of the cube, essentially reducing to two dimensional assembly. Orientations do not matter in 2D assembly, hence a default value is used. 

## On splitting across parallel runs:
`tot_splits` split the total number of input rules into separate batches and stores the valid outputs and their frequencies, complexities separately (referred to as `nsplit`). Later, you can use `combine.py` or similar file to identify unique shapes and add up their frequencies and store as combined_files. 
