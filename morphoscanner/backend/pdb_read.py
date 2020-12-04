"""
Created on Wed Nov 18 19:55:01 2020

@author: lillo
"""
import numpy as np
from timeit import default_timer as timer
import backend


from .distance import get_euclidean_distance
from .distance_tensor import get_coordinate_tensor_from_dict_single, fast_cdist


def get_coordinate_from_pdb(file):
    '''
    Parse a pdb file. Support single chain and multiple chain

    Parameters
    ----------
    file : str
        The path of the .pdb file in your system.

    Returns
    -------
    coordinate_dict : dict
        A dict of dict with the coordinate of each atom of the pdb file.
        
        Depending on the input file it has different levels of nesting:
            
            for single chain:
                atom_index : [x,y,z]
                
            for multiple chain:
                
                chain_index : {atom index : [x,y,z]}
    '''
    with open(file) as pdbfile:

        coordinate_dict = {}
        atom_count_dict = {}
        start = 0

        for line in pdbfile:
            
            # split line
            splitted_line = [line[:6], line[6:11], line[12:16], line[17:20], line[21], line[22:26], line[30:38], line[38:46], line[46:54]]
            # get line header
            line_id = splitted_line[0].split()[0]
            
            #check for atom and heteroatom
            if line_id in {'ATOM', 'HETATM'}:
                
                # get CA atom only
                if splitted_line[2].split()[0] in {'CA'}:
                    
                    # get atom num for indexing
                    #atom_num = int(splitted_line[5])
                    # get protein chain for indexing
                    chain = splitted_line[4]
                    # get coordinates
                    x, y, z = float(splitted_line[6]), float(splitted_line[7]), float(splitted_line[8])
                    
                    # check if actual chain already has an entry in coordinate_dict
                    if chain not in coordinate_dict.keys():
                        
                        # index from 'start'
                        atom_count_dict[chain] = start
                        # create key for new chain
                        coordinate_dict[chain] = {}
                        # put actual atom coordinates in coordinate_dict
                        coordinate_dict[chain][atom_count_dict[chain]] = np.array([x,y,z])
                    # if actual chain already in coordinate_dict
                    else:
                        # move index forward
                        atom_count_dict[chain] += 1
                        # add the atom coordinates
                        coordinate_dict[chain][atom_count_dict[chain]] = np.array([x,y,z])

    # if there is only one chain, flat the dict
    if len(coordinate_dict) == 1:
        coordinate_dict = coordinate_dict.get([k for k in coordinate_dict][0])

    return coordinate_dict


def compute_distance_map_protein(coordinate_dict):
    '''Compute pointwise distance map on a dict of coordinate

    Parameters
    ----------
    coordinate_dict : dict
        The dict of coordinate, in the form {atom : [x, y, z]}.

    Returns
    -------
    distance_map : numpy.ndarray
        numpy.ndarray of shape (len(coordinate_dict), len(coordinate_dict))
        each element i,j is the distance between the corresponding elements in
        coordinate_dict[i] and coordinate_dict[j]
    '''
    # create zeros map
    distance_map = np.zeros((len(coordinate_dict),len(coordinate_dict)))
    # iterate through each dict element
    for value1 in coordinate_dict:
        # get element coordinate
        coordinate_1 = coordinate_dict[value1]
        # iterate in the upper triangle
        # skip main diagonal because is 0 already (np.zeros)
        for value2 in range(value1+1, len(coordinate_dict)):
            # get secod element coordinate
            coordinate_2 = coordinate_dict[value2]
            # calculate distance
            euclidean_distance = get_euclidean_distance(coordinate_1, coordinate_2)
            # fill array with distance
            distance_map[value1][value2] = euclidean_distance
            # fill lower triangle with distance
            distance_map[value2][value1] = euclidean_distance
            
    return distance_map


def test_distance_error(protein_dict):
    
    # print number of CA
    print('Your protein has %d CA atoms.\n' % len(protein_dict))
    
    # compute pointwise distance
    start_correct = timer()
    correct_distance = compute_distance_map_protein(protein_dict)
    end_correct = timer()
    correct_time = (end_correct-start_correct)
    print('Correct calculation of distances needs %s seconds.\n' % str(correct_time))
    
    # parallel distance
    start_init_tesor = timer()
    protein_tensor = get_coordinate_tensor_from_dict_single(protein_dict)
    end_init_tensor = timer()
    init_tens_time = (end_init_tensor-start_init_tesor)
    print('Time to initialize tensor was %s seconds.' % str(init_tens_time))
    
    #c compute distance by matrix operations
    start_tens_multi = timer()
    approximate_distance = fast_cdist(protein_tensor, protein_tensor)
    end_tens_multi = timer()
    tens_multi_time = (end_tens_multi-start_tens_multi)
    print('Time to compute tensor operation was %s seconds.' % str(tens_multi_time))
    print('Total time to compute tensor operation was %s seconds.\n' % str(init_tens_time+tens_multi_time))
    #calculate error matrix
    error = correct_distance - approximate_distance.cpu().numpy()
    
    # print number of distances on which statistics are calculated
    print('Statistics are calculated on %d measurements.\n' % error.size)
    
    # calculate and print statistics
    mae = np.mean(np.abs(error))
    print('Mean absolute error (MAE) is:', mae)
    max_error = np.max(error)
    print('Maximum error is:', max_error)
    median_error = np.median(error)
    print('Median error is:', median_error)
    mse = np.mean(np.square(error))
    print('Mean Squared _Error is:', mse)
    std = np.std(error)
    print('Standard deviation is: ', std)
    
    return


def get_distance_map_from_pdb(path:str):
    coordinate = backend.pdb_read.get_coordinate_from_pdb(path)
    keys = [k for k in coordinate.keys()]
    print('The file you opened has the following models ID:\n')
    print(keys)
    model = input('\nWrite the models that you want to select:')
    if len(model)>0:
        try:
            model_coord = coordinate[model]
            print('\nModel %s was selected' % str(model))
            tensor = backend.pdb_read.get_coordinate_tensor_from_dict_single(model_coord)
            distance_map = backend.distance_tensor.fast_cdist(tensor,tensor)
        except KeyError:
            print('%s is not a correct model ID.' % model)
    else:
        model = keys[0]
        print('\nModel %s was selected.' % str(model))
        tensor = backend.pdb_read.get_coordinate_tensor_from_dict_single(coordinate[model])
        distance_map = backend.distance_tensor.fast_cdist(tensor,tensor)
    return distance_map