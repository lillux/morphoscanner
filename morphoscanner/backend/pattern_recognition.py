"""
Functions for pattern recognition of beta sheet

    class cross_correlation():

@author: lillo
"""

from .utilities import contact_list_from_dict

import numpy as np
import pandas as pd
import tqdm


# COMPUTE CONTACT MAPS
# TO DO: parametrize the threshold distance in a better way (e.g. )
# DISTANCE MAP ARRAY WILL NOT EXIST IN DISTINCT PEPTIDE SIZE SIMULATION
def compute_contact_maps_as_array(distance_maps_array, radius_multiplier=1.5):      ## TAKE INPUT FOR RADIUS MULTIPLIER

    # median consecutive residues distance from first peptide
    distances_pep_1 = []
    for i in range(distance_maps_array[0][0].shape[0] - 1): # there was a -1, and len(denoised_dict) was len(coordinate_dict)-1
        distances_pep_1.append(distance_maps_array[0][0][i][i+1])

    intrapeptide_minimum_distance = np.median(distances_pep_1)

    #temporary list
    contact_map_list = []

    # contact is in a distance up to 150% of the intrapeptide_minimum_distance [TO IMPROVE!!!]
    threshold_distance = (intrapeptide_minimum_distance * radius_multiplier)

    for model_1 in range(distance_maps_array.shape[0]):
        contact_map_list.append([])
        for model_2 in range(distance_maps_array[model_1].shape[0]):

            contact_map_list[model_1].append([])

            if model_1 == model_2:

                contact_map_list[model_1][model_2].extend(np.zeros((distance_maps_array.shape[2], distance_maps_array.shape[3])))

            else:

                contact_map = np.zeros((distance_maps_array.shape[2], distance_maps_array.shape[3]))

                for chain_1 in range(distance_maps_array[model_1][model_2].shape[0]):

                    for chain_2 in range(distance_maps_array[model_1][model_2][chain_1].shape[0]):

                        distance = distance_maps_array[model_1][model_2][chain_1][chain_2]

                        if distance < threshold_distance:
                            contact_map[chain_1][chain_2] = 1 #True
                        else:
                            pass

                contact_map_list[model_1][model_2].extend(contact_map)

    contact_array = np.asarray(contact_map_list)

    return contact_array

 

        
#### ANALYSIS

def shift_library_maker(contact_map_to_analyze):

    ''' Create shift matrix library to perform pattern recognition on a contact map.
    
    Inputs: numpy.array, a single contact map
    
    returns: dict of shift matrix to analyze the given contact map

            shift matrix are diagonal matrix spanning all the given contact map
            
    '''
    
    row = contact_map_to_analyze.shape[0]
    col = contact_map_to_analyze.shape[1]

    kron_dict = {}
    kron_list_parallel = []
    kron_list_antiparallel = []

    for e in range(-row+1, col):
        array = np.eye(row, col, e)
        kron_list_parallel.append(array)
        kron_list_antiparallel.append(np.fliplr(array))

    kron_array_parallel = np.asarray(kron_list_parallel)
    kron_array_antiparallel = np.asarray(kron_list_antiparallel)

    kron_dict['parallel'] = kron_array_parallel
    kron_dict['antiparallel'] = kron_array_antiparallel

    return kron_dict




def normalized_cross_correlation_function(contact_map, minimum_contact=2):
    '''
    Calculate normalized cross correlation function between a contact map and an ideal map.

    Arguments : contact map, as output from get_contact_maps function
                shift_matrix_stack, as output from shift_matrix_maker function

    Return : a list [ncc_value, index (in the shift_matrix_stack) of the shift matrix
                that is matching the contact map

            '''
    shift_matrix_library = shift_library_maker(contact_map)

    cross_correlation_values = []
    max_val = []
    sum_contact_map = np.sum(contact_map)

    if sum_contact_map < minimum_contact:
        pass

    else:
        for sense in shift_matrix_library:
            for index, z in enumerate(shift_matrix_library[sense]):

                shift_matrix = shift_matrix_library[sense][index]
                sum_shift_matrix = np.sum(shift_matrix)
                ncc_value = (np.sum((contact_map * shift_matrix))/((np.sqrt(sum_contact_map))*(np.sqrt(sum_shift_matrix))))  # normalized cross correlation function of contact matrix and shift matrix
                cross_correlation_values.append([ncc_value, index, sum_contact_map, sense])

            max_val = max(cross_correlation_values) # get only the best match (highest value of ncc)

    return max_val




def normalized_cross_correlation_for_dataset(contact_array):
    '''Calculate normalized cross correlation function between the full contacts map and
    the shift matrix.

    Arguments : contact map, as output from get_contact_maps function
                shift_matrix_stack, as output from shift_matrix_maker function

    Return : a list [ncc_value, index (in the shift_matrix_stack) of the shift matrix that is matching the contact map'''

    contact_dict = {}

    #for row in tqdm.tqdm(range(contact_array.shape[0])):
    for row in range(contact_array.shape[0]):

        for col in range((row+1), contact_array.shape[1]):
        #for col in range(contact_array.shape[1]):

            best_match = []
            best_match = normalized_cross_correlation_function(contact_array[row][col])

            if len(best_match) == 0:
                pass

            else:
                if row in contact_dict:
                    contact_dict[row].append([row, col, best_match])

                else:
                    contact_dict[row] = [[row, col, best_match]]

    return contact_dict



#contact_array = frame_contact
def cross_correlation_function_for_dataset_with_dataframe(contact_array):
    '''Perform Normalized Cross Correlation function on the dataset
        to check for contact. Get a dict for processing and a pandas.DataFrame
        for data analysis

        Input: contact maps

        Output: contact_dict,         for further processing
                pandas.DataFrame,     for data analysis


    '''
    contact_dict = {}

    for row in range(contact_array.shape[0]):

        for col in range((row+1), contact_array.shape[1]):
            best_match = []
            best_match = normalized_cross_correlation_function(contact_array[row][col])

            if len(best_match) == 0:
                pass

            else:
                if row in contact_dict:
                    contact_dict[row].append([row, col, best_match])

                else:
                    contact_dict[row] = [[row, col, best_match]]

    contact_list = contact_list_from_dict(contact_dict)

    columns_names = ['peptide1', 'peptide2', 'NCC Value', 'shift index', 'contacts', 'sense']

    df = pd.DataFrame(contact_list, columns=columns_names)

    return contact_dict, df




#denoise dataset. GET SENSE 'PARALLEL' OR 'ANTIPARALLEL'....'NEED TO KNOW AMINOACID OF PEPTIDES'
def denoise_contact_maps(contact_maps):
    
    '''Denoise the contact_maps dataset using the shift_matrix
    
    Arguments : contact_maps, normalized_cross_correlation_result
    
    return : a dict with key:value = row : row, col, denoised_map
    
    '''

    normalized_cross_correlation_results, df = cross_correlation_function_for_dataset_with_dataframe(contact_maps)


    denoised_dict = {}

    for peptide_1 in normalized_cross_correlation_results:
        denoised_dict[peptide_1] = {}
        for index, peptide_2 in enumerate(normalized_cross_correlation_results[peptide_1]):

            row = peptide_2[0]
            col = peptide_2[1]



            contact_map = contact_maps[row][col]
            sense = peptide_2[2][3]
            shift_matrix_index = normalized_cross_correlation_results[peptide_1][index][2][1]

            shift_matrix = shift_library_maker(contact_map)
            shift_matrix = shift_matrix[sense][shift_matrix_index]
            denoised_map = contact_map * shift_matrix

            denoised_dict[row][col] = denoised_map
            
            
    full_denoised_dict = {}
    for peptide_1 in tqdm.tqdm(denoised_dict):
        for peptide_2 in denoised_dict[peptide_1]:
            contact_map = denoised_dict[peptide_1][peptide_2]

            if peptide_1 in full_denoised_dict:
                full_denoised_dict[peptide_1][peptide_2] = contact_map

            if peptide_1 not in full_denoised_dict:
                full_denoised_dict[peptide_1] = {peptide_2:contact_map}

            if peptide_2 in full_denoised_dict:
                full_denoised_dict[peptide_2][peptide_1] = contact_map.T

            if peptide_2 not in full_denoised_dict:
                full_denoised_dict[peptide_2] = {peptide_1:contact_map.T}
    
    return full_denoised_dict, df

