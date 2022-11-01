"""
    Functions for common functionality of backend


@author: lillo
"""
import torch

def get_coordinate_dict_from_array(array):
    '''Take decomposition result and convert it into a coordinate vectors dict

    Argument: decomposition results

    return: dict with reconstructed 3d coordinate vector

    '''

    # make list with n == array.shape[0]
    coordinate = [e for e in array]

    # initialize empty dict
    reconstructed_coordinate_dict = {}

    # fill the dict with the ccordinate vectors
    for index,coordinate_vector in enumerate(coordinate):
        reconstructed_coordinate_dict[index] = coordinate_vector

    return reconstructed_coordinate_dict


#fixed for denoised
def contact_list_from_dict_v1(contact_dict):
    
    contact_list = []
    for peptide in contact_dict:

        for contact in contact_dict[peptide]:

            new_data = [contact[0], contact[1], contact[2][0], contact[2][1], contact[2][2], contact[2][3], contact[2][5]]
            contact_list.append(new_data)
    return contact_list


# return string in a string with numbers
def return_if_string(string):
    digits = []
    for i in string:
        if not i.isdigit():
            digits.append(i)

    string = ''.join(digits)

    return string


# return numbers in a string with numbers
def return_if_digit(string):
    digits = []
    for i in string:
        if i.isdigit():
            digits.append(i)

    string = ''.join(digits)

    return string

def get_map_index(self, frame=0):
    '''
    Compute index for fast data retrival from a single distance map, or tensor, using indexing
    '''
    # start from 0
    total = 0
    # instantiate empty dict
    index = {}
    # for each parsed peptide
    for peptide in self.frames[frame].peptides:
        # calculate number of peptide grains
        len_pep = len(self.frames[frame].peptides[peptide].coordinates)
        # save [start, end] index values
        index[peptide] = [total,(total+len_pep)]
        # add actual peptide number of grains to total
        total += len_pep
    # save as instance attibute
    # self.map_index = index # use it to make instance attribute
    return index

def retrieve_map_from_distance_map(dist_map: torch.tensor, map_index: dict, i: int, j: int):
    ij_map = dist_map[map_index[i][0]:map_index[i][1], map_index[j][0]:map_index[j][1]]
    return ij_map