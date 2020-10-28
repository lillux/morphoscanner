"""
Functions for pattern recognition of beta sheet

    class cross_correlation():

@author: lillo
"""

from .utilities import contact_list_from_dict_v1

import torch
import numpy as np
import pandas as pd
import tqdm

### Torch functions

def shift_library_maker_torch(contact_map_to_analyze, device = 'cpu'):
    ''' 
    Create shift matrix library to perform pattern recognition on a contact map.
    Shift matrix are diagonal matrix spanning all the given contact map. 

    Arguments : numpy.array, a single contact map of size n*m
    
    returns : dict of torch.tensor
            
        These are the shift matrix to analyze the given contact map

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

    kron_dict['parallel'] = torch.from_numpy(np.asarray(kron_list_parallel)).to(device)
    kron_dict['antiparallel'] = torch.from_numpy(np.asarray(kron_list_antiparallel)).to(device)
    return kron_dict
        
        
## WORKING GOOD IN PARALLEL
def normalized_cross_correlation_function_torch_v1(contact_map, minimum_contact=2, device='cpu'):
    '''
    Calculate normalized cross correlation function between a contact map and an ideal map.

    Arguments : contact map, as output from get_contact_maps function
                shift_matrix_stack, as output from shift_matrix_maker function

    return : a list [ncc_value, index (in the shift_matrix_stack) of the shift matrix
                that is matching the contact map
    '''
    contact_map = contact_map.double()
    shift_matrix_library = shift_library_maker_torch(contact_map, device=device)

    cross_correlation_values = []
    max_val = []
    sum_contact_map = torch.sum(contact_map)
    shift_matrix_center_index = ((contact_map.shape[0] + contact_map.shape[1]) -1)//2

    if sum_contact_map < minimum_contact:
        pass

    else:
        
        for sense in shift_matrix_library:
            
            signal_full = contact_map * shift_matrix_library[sense]
            signal_tens = torch.sum(signal_full, dim=(1,2))
            norm = torch.sqrt(sum_contact_map) * torch.sqrt(torch.sum(shift_matrix_library[sense], dim=(1,2)))
            ncc = signal_tens/norm
            ncc_index = torch.argmax(ncc)
            ncc_val = ncc[ncc_index]
        
            denoised = signal_full[ncc_index]
            cross_correlation_values.append([ncc_val, ncc_index, sum_contact_map, sense, denoised])
            
        max_val = max(cross_correlation_values) # get only the best match (highest value of ncc)
        sum_denoised = torch.sum(max_val[4])
        if sum_denoised >= minimum_contact:
            
            shift = shift_matrix_center_index - max_val[1]
            max_val[2] = sum_denoised
            max_val.append(shift)
            
            return max_val
        
        else:
            pass


def cross_correlation_function_for_dataset_with_dataframe_torch_v1(contact_array, device='cpu'):
    '''
    Perform Normalized Cross Correlation function on the dataset
    to check for contact. Get a dict for processing and a pandas.DataFrame
    for data analysis

        Input: contact maps

        Output: contact_dict,         for further processing
                pandas.DataFrame,     for data analysis
    '''
    contact_dict = {}
    denoised_dict = {}

    for row in range(len(contact_array)):

        for col in range((row+1), len(contact_array[row])):
            best_match = []

            best_match = normalized_cross_correlation_function_torch_v1(contact_array[row][col], device=device)

            if best_match == None:
                pass

            else:
                if row in contact_dict:
                    contact_dict[row].append([row, col, best_match])

                else:
                    contact_dict[row] = [[row, col, best_match]]

                if row in denoised_dict:
                    denoised_dict[row][col] = best_match[4]

                else:
                    denoised_dict[row] = {col : best_match[4]}

    contact_list = contact_list_from_dict_v1(contact_dict)
    columns_names = ['peptide1', 'peptide2', 'NCC Value', 'shift index', 'contacts', 'sense', 'shift']

    df = pd.DataFrame(contact_list, columns=columns_names)

    return contact_dict, df, denoised_dict


def denoise_contact_maps_torch_v1(contact_maps, device='cpu'):
    
    '''Denoise the contact_maps dataset using the shift_matrix
    
    Arguments : contact_maps, normalized_cross_correlation_result
    
    return : a dict with key:value = row : row, col, denoised_map
    
    '''
    normalized_cross_correlation_results, df, denoised_dict = cross_correlation_function_for_dataset_with_dataframe_torch_v1(contact_maps, device=device)

    full_denoised_dict = {}
    for peptide_1 in tqdm.tqdm(denoised_dict):

        for peptide_2 in denoised_dict[peptide_1]:
            contact_map = denoised_dict[peptide_1][peptide_2]

            if peptide_1 in full_denoised_dict:
                full_denoised_dict[peptide_1][peptide_2] = contact_map

            if peptide_1 not in full_denoised_dict:
                full_denoised_dict[peptide_1] = {peptide_2:contact_map}

            if peptide_2 in full_denoised_dict:
                #full_denoised_dict[peptide_2][peptide_1] = contact_map.T #this is for numpy
                full_denoised_dict[peptide_2][peptide_1] = contact_map.transpose(1,0) # this is for pytorch


            if peptide_2 not in full_denoised_dict:
                #full_denoised_dict[peptide_2] = {peptide_1:contact_map.T} #this is for numpy
                full_denoised_dict[peptide_2] = {peptide_1:contact_map.transpose(1,0)} # this is for pytorch
    
    return full_denoised_dict, df