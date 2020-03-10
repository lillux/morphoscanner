# -*- coding: utf-8 -*-


import numpy as np

import re



class distance():
        
        
        # compute euclidean distance
        def get_euclidean_distance(point_1, point_2):

            euclidean_distance = np.sqrt(np.sum([((point_1[0] - point_2[0])**2), ((point_1[1] - point_2[1])**2), ((point_1[2] - point_2[2])**2)]))

            return euclidean_distance
        
        

        # compute distance map between two peptides
        def compute_distance_map(coordinate_dict, peptide_1, peptide_2):

            distance_map = []
            for amino_1 in coordinate_dict[peptide_1]:
                coordinate_1 = coordinate_dict[peptide_1][amino_1]

                distance_map.append([amino_1])

                for amino_2 in coordinate_dict[peptide_2]:
                    coordinate_2 = coordinate_dict[peptide_2][amino_2]

                    euclidean_distance = morphoscanner.distance.get_euclidean_distance(coordinate_1, coordinate_2)
                    distance_map[amino_1].append(euclidean_distance)

                del distance_map[amino_1][0]

            distance_map = np.asarray(distance_map)

            return distance_map

        
        
        # compute distance map and return a n_peptide x n_peptide x n_res x n_res array
        def compute_distance_maps_from_coordinate_dict(coordinate_dict):

            aggregate_distance_map = []

            for peptide_1 in tqdm.tqdm(coordinate_dict):
            #for peptide_1 in coordinate_dict:
                aggregate_distance_map.append([peptide_1])

                #for peptide_2 in tqdm.tqdm(coordinate_dict):
                for peptide_2 in coordinate_dict:
                    distance_map = morphoscanner.distance.compute_distance_map(coordinate_dict, peptide_1, peptide_2)

                    aggregate_distance_map[peptide_1].append(distance_map)

                del aggregate_distance_map[peptide_1][0]

            aggregate_distance_array = np.asarray(aggregate_distance_map)

            return aggregate_distance_array


        # COMPUTE CONTACT MAPS
        # TO DO: parametrize the threshold distance in a better way (e.g. )
        def compute_contact_maps_as_array(distance_maps_array):

            # distance between the first and the second aminoacid of the first chain
            intrapeptide_minimum_distance = distance_maps_array[0][0][0][1] 

            contact_map_list = []

            # contact is in a distance up to 150% of the intrapeptide_minimum_distance [TO IMPROVE!!!]
            threshold_distance = (intrapeptide_minimum_distance * 1.5)

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
        
        
        
        # get average distance map from distance maps set
        def get_mean_distance_map(distance_maps):
            '''
            Calculate mean distance map from distance maps set

            Argument: distance maps set

            return: np.array with average intrapeptide distance

            '''

            # create array of zeros of shape number_of_residues * number_of_residues
            # depending on peptide residue number ### TO FIX FOR MULTIMONOMERIC ASSEMBLY
            base = np.zeros((distance_maps[0][0].shape[0], distance_maps[0][0].shape[1]))

            # initialize counter
            counter = 0

            # iterate throught peptides in the aggregate
            for peptide_1 in range(distance_maps.shape[0]):
                for peptide_2 in range(distance_maps.shape[1]):

                    # if peptide index are the same (intrapeptide distance map)
                    if peptide_1 == peptide_2:

                        # intrapeptide distance map
                        actual_distance_map = distance_maps[peptide_1][peptide_2]

                        # sum base and current distance map
                        base = base + actual_distance_map

                        #update counter
                        counter += 1

            #for element in base (every element is the sum of distance_map(i,j) for every distance map)
            for row in range(len(base)):
                for col in range(len(base)):

                    # find the mean for every element of the cumulative distance map
                    base[row][col] = (base[row][col])/counter

            return base

        # Singolar Value Decomposition of distance_map
        def decompose_distance_map(distance_map):
            '''Use Singular value decomposition to get

            distance_map.shape[1] dimensional coordinate
            (same n of dimension as the peptide n of residue)

            As described in:
            Mathematical Modeling of Protein Structure Using Distance Geometry
            Jeong-Mi Yoon, Yash Gad, Zhijun Wu

            Argument: distance map (numpy.array 2D)
            return: X : actual decomposition


            '''

            # initialize a zeros matrix of same shape as the input map
            D = np.zeros(distance_map.shape)

            #iterate trought row
            for i in range(distance_map.shape[0]):

                # iterate trought cols
                for j in range(distance_map.shape[1]):

                    # distance between point point i and point j 
                    dij = distance_map[i][j]

                    # distance between point 0 and point j
                    d0j = distance_map[0][j]

                    #distance between point i and point 0
                    di0 = distance_map[i][0]

                    #fill the zeros matrix with the value obtained with this formula
                    D[i][j] = (d0j**2 + di0**2 - dij**2)/2

            # check rank of matrix (should be of rank 3, but it is of rank distance_map.shape[1])
            #rank = np.linalg.matrix_rank(D)

            # Singular value decomposition on the D matrix
            #svd = np.linalg.svd(D)

            svd = np.linalg.svd(D, full_matrices=False)

            # Calculate distance_map.shape[1] dimensional coordinate, but you need 3
            # the non necessary dimension can give data to better reconstruct the peptide structure
            X = svd[0]*np.sqrt(svd[1])


            #return X, svd, D, rank
            return X

        def get_coordinate_from_decomposition(decomposition):
            '''Take decomposition result and convert it into a coordinate vectors dict

            Argument: decomposition results

            return: dict with reconstructed 3d coordinate vector

            '''

            # take only the first three value to compose a 3D coordinate vector
            coordinate = [e[:3] for e in decomposition]

            # initialize empty dict
            reconstructed_coordinate_dict = {}

            # fill the dict with the ccordinate vectors
            for index,coordinate_vector in enumerate(coordinate):
                reconstructed_coordinate_dict[index] = coordinate_vector

            return reconstructed_coordinate_dict


        # reconstruct 3d coordinate from a distance map
        def get_coordinate_from_distance_map(distance_map):
            ''' compute 3d coordinate from distance map

            Argument: distance_map (numpy.array)

            return: dict with 3d coordinate for every alpha-carbon of a peptide

            '''
            # perform singular value decomposition on distance_map (preprocessed)
            decomposed_mean_distance_map = decompose_distance_map(distance_map)


            # get 3D coordinate
            reconstructed_coordinate_dict = get_coordinate_from_decomposition(decomposed_mean_distance_map)

            return reconstructed_coordinate_dict
        








































