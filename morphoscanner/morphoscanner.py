#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 01:01:39 2020

@author: lillo
"""


import math
import itertools
import numpy as np
import matplotlib.pyplot as plt

import plotly
import plotly.express as px

from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
import tqdm
#from functools import lru_cache
#import re
import networkx as nx
from networkx.algorithms import approximation


import MDAnalysis as mda

#import scipy
#import sklearn
#import skimage

#import xml.etree.ElementTree as et
#from Bio.PDB import *
import nglview as nv

from timeit import default_timer as timer



# http://nglviewer.org/nglview/latest/api.html
# https://biopython.org/wiki/The_Biopython_Structural_Bioinformatics_FAQ
# https://ambermd.org/tutorials/analysis/tutorial_notebooks/nglview_notebook/index.html
# https://amber-md.github.io/pytraj/latest/_api/pytraj.html

#####################################################################
###############################################################################
##################################################################

#####################  IMPORT

#contact_matrix = np.loadtxt('/home/lillo/TesiCNTE/CNTE/dataset/contact_matrix.txt')   #laptop
#contact_matrix = np.loadtxt('/home/lillo/Code/Tesi/dataset/contact_matrix.txt')        #fisso
#contact_matrix_single = contact_matrix.reshape(100,100,12,12)

#gromacs_output = open('/home/lillo/Code/Tesi/dataset/dm4500Compl_mix1_K2_1%4500ns.gro') #fisso
#gromacs_output = open('/home/lillo/TesiCNTE/CNTE/dataset/dm4500Compl_mix1_K2_1%4500ns.gro') #laptop

path = '/home/lillo/Code/Tesi/dataset/dm4500Compl_mix1_K2_1%4500ns.gro' #fisso
#path = '/home/lillo/TesiCNTE/CNTE/dataset/dm4500Compl_mix1_K2_1%4500ns.gro' #laptop

# import 2mxu file (beta sheet)

#path_to_mmCIF = open('/home/lillo/TesiCNTE/pdb/2mxu/2mxu.cif')  ## laptop
#path_to_pdb = '/home/lillo/TesiCNTE/pdb/2mxu/2mxu.pdb'  ## laptop
#pa_to_pdb = '/home/lillo/TesiCNTE/pdb/2mxu/2mxu.pdb'  ## laptop

#path_to_mmCIF = open('/home/lillo/Code/Tesi/pdb/2mxu/2mxu.cif')  ## fisso
#path_to_pdb = '/home/lillo/Code/Tesi/pdb/2mxu/2mxu.pdb'  ## fisso
#pa_to_pdb = '/home/lillo/Code/Tesi/pdb/2mxu/2mxu.pdb'  ## fisso

# aggregate blob

#seed_1_path = '/home/lillo/TesiCNTE/from_cluster/aggregate1.gro' # laptop
seed_1_path = '/home/lillo/Code/Tesi/dataset/aggregate1.gro' # Fisso

# Trajectory with aggregate seed
#trj_xtc = '/home/lillo/TesiCNTE/CNTE/trajectory/prd-LDLK12-100mer-out-mol.xtc'  #laptop
#trj_gro = '/home/lillo/TesiCNTE/CNTE/trajectory/min-LDLK12-100mer-out-c.gro'  #laptop

trj_xtc = '/home/lillo/Code/Tesi/dataset/trajectory_6_12_19/prd-LDLK12-100mer-out-mol.xtc'  #fisso
trj_gro = '/home/lillo/Code/Tesi/dataset/trajectory_6_12_19/min-LDLK12-100mer-out-c.gro'  #fisso



########################################################
##################################################################################
##################################################################

class morphoscanner():
    
    
    
    
    
    
    
    
    class read():
        
        # READ .gro FILE AND PREPROCESSING

        def clean_gro(path):


                # open file .gro and return a list with one element per line of the .gro file
            def read_gro(path):
                gromacs_output = open(path)

                gro_file = []
                for line in tqdm.tqdm(gromacs_output):
                    gro_file.append(line)



                gromacs_output.close()        

                return gro_file



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


            # remove first, second and last lines from gro_file and reorder information
            # FIX OPTION TO GET ENTRY RELATED TO A LABEL (as 'bb' or 'ca')
            def clean_gro_file(gro_file):
                cleaned_gro_file = []
                for aminoacid in tqdm.tqdm(gro_file[2:-1]):
                    splitted = aminoacid.split()
                    if splitted[1] == 'BB':
                        position_in_peptide = return_if_digit(splitted[0])
                        residue = return_if_string(splitted[0])
                        index = splitted[2]
                        x = splitted[3]
                        y = splitted[4]
                        z = splitted[5]
                        cleaned_gro_file.append([index, position_in_peptide, residue, x, y, z])
                return cleaned_gro_file


            gro_file = read_gro(path)
            cleaned_gro_file = clean_gro_file(gro_file)

            return cleaned_gro_file
        
        # create coordinate dict from cleaned_gro_file
        def get_coordinate_dict_from_cleaned_gro(cleaned_gro_file):

            peptide_lenght_list = []

            temporary_list = []

            # iterate trough cleaned_gro_file
            for residue in cleaned_gro_file:

                # if temporary list just started, add aminoacid position in chain
                if len(temporary_list) == 0:
                    temporary_list.append(int(residue[1]))

                else:
                    # if position of actual residue is less than last residue
                    if temporary_list[-1] > int(residue[1]):

                        # append lenght of last peptide to peptide lenght list
                        peptide_lenght_list.append(len(temporary_list))

                        # empty temporary list
                        temporary_list = []

                        # append actual residue position
                        temporary_list.append(int(residue[1]))

                    # if position of actual residue is higher than last residue, ad current residue position
                    else:
                        temporary_list.append(int(residue[1]))

            # append last peptide lenght to lenght stack
            peptide_lenght_list.append(len(temporary_list))

            # create empty dict for coordinate
            peptide_coordinate_dict = {}

            # create an entry in dict for every peptide in the file
            for peptide in range(len(peptide_lenght_list)):
                peptide_coordinate_dict[peptide] = {}

                # for every residue in lenght peptide, add coordinate x, y, z
                for residue in range(peptide_lenght_list[peptide]):
                    peptide_coordinate_dict[peptide][residue] = [float(coordinate) for coordinate in cleaned_gro_file[(peptide * peptide_lenght_list[peptide])+residue][3:]]

            return peptide_coordinate_dict
        
        
        def get_coordinate_dict_from_cleaned_gro_for_fixed_lenght_peptides(cleaned_gro_file, peptide_lenght):
            '''Works only with system made of peptide of equal peptide length
                (the peptide in the trajectory have all the same number of aminoacids).

            Usefull to disassemble premade aggregate in a system of omogenous peptide lenght'''

            peptide_coordinate_dict = {}

            for peptide in range(len(cleaned_gro_file)//peptide_lenght):

                peptide_coordinate_dict[peptide] = {}

                peptide_data = cleaned_gro_file[(peptide*peptide_lenght):((peptide*peptide_lenght)+peptide_lenght)]

                for index, amino in enumerate(peptide_data):

                    peptide_coordinate_dict[peptide][index] = [float(amino[3]), float(amino[4]), float(amino[5])] ## can you generalize this 3? Maybe with re or isdigit function??

            return peptide_coordinate_dict
        

    
    
    
    
    
    
    
    
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
        
        
        
        
        
        
        
        
        
    class cross_correlation():
        
        #### ANALYSIS

        def shift_library_maker(contact_map_to_analyze):

            ''' riceve numero di righe e di colonne
            restituisce un array shape((((row + col)*2)-2),row,col).
            ogni slice è una diagonale. Lo stack copre le diagonali su tutta la matrice'''

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


        def normalized_cross_correlation_function(contact_map):
            '''
            Calculate normalized cross correlation function between a contact map and an ideal map.

            Arguments : contact map, as output from get_contact_maps function
                        shift_matrix_stack, as output from shift_matrix_maker function

            Return : a list [ncc_value, index (in the shift_matrix_stack) of the shift matrix
                        that is matching the contact map

                    '''
            shift_matrix_library = morphoscanner.cross_correlation.shift_library_maker(contact_map)

            cross_correlation_values = []
            max_val = []
            sum_contact_map = np.sum(contact_map)

            if sum_contact_map < 2:
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
            '''Calculate normalized cross correlation function between the full contacts map and and the .

            Arguments : contact map, as output from get_contact_maps function
                        shift_matrix_stack, as output from shift_matrix_maker function

            Return : a list [ncc_value, index (in the shift_matrix_stack) of the shift matrix that is matching the contact map'''

            contact_dict = {}

            #for row in tqdm.tqdm(range(contact_array.shape[0])):
            for row in range(contact_array.shape[0]):

                for col in range((row+1), contact_array.shape[1]):
                #for col in range(contact_array.shape[1]):

                    best_match = []
                    best_match = morphoscanner.cross_correlation.normalized_cross_correlation_function(contact_array[row][col])

                    if len(best_match) == 0:
                        pass

                    else:
                        if row in contact_dict:
                            contact_dict[row].append([row, col, best_match])

                        else:
                            contact_dict[row] = [[row, col, best_match]]

            return contact_dict
        
        
        
        
        
    
    class denoise():
        
        
        #denoise dataset
        def denoise_contact_maps(contact_maps):
            
            '''Denoise the contact_maps dataset using the shift_matrix
            
            Arguments : contact_maps, normalized_cross_correlation_result
            
            return : a dict with key:value = row : row, col, denoised_map
            
            '''
        
            normalized_cross_correlation_results = morphoscanner.cross_correlation.normalized_cross_correlation_for_dataset(contact_maps)
        
        
            denoised_dict = {}
        
            for peptide_1 in normalized_cross_correlation_results:
                denoised_dict[peptide_1] = {}
                for index, peptide_2 in enumerate(normalized_cross_correlation_results[peptide_1]):
        
                    row = peptide_2[0]
                    col = peptide_2[1]
        
        
        
                    contact_map = contact_maps[row][col]
                    sense = peptide_2[2][3]
                    shift_matrix_index = normalized_cross_correlation_results[peptide_1][index][2][1]
        
                    shift_matrix = morphoscanner.cross_correlation.shift_library_maker(contact_map)
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
            
            
            return full_denoised_dict
        
        
        












