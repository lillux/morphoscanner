"""
    Functions for common functionality of backend


@author: lillo
"""



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

        
        
def contact_list_from_dict(contact_dict):
    
    contact_list = []
    for peptide in contact_dict:

        for contact in contact_dict[peptide]:

            new_data = [contact[0], contact[1], contact[2][0], contact[2][1], contact[2][2], contact[2][3]]
            contact_list.append(new_data)
    return contact_list