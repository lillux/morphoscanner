"""
Created on Wed Nov 18 17:50:56 2020

@author: lillo
"""
import numpy as np


def contact_map_helix_torch(distance_map):
    
    t_h_sub = distance_map.clone()
    t_h_sub[t_h_sub <= 4.4] = 0
    t_h_sub[(t_h_sub > 4.4) & (t_h_sub<=5.1)] = 1
    t_h_sub[(t_h_sub > 5.1) & (t_h_sub<=6.5)] = 2
    t_h_sub[t_h_sub > 6.5] = 3
    
    return t_h_sub


def contact_tracer(contact_map):
    contact_map = contact_map.clone()
    pace_helix = 0
    count_alpha = 0
    percentage_alpha = 0
    count_beta  = 0
    percentage_beta = 0
    length_beta = 0
    folding = {}
    explored = []
    for i in range(contact_map.shape[0]):
        for j in range(i+1, contact_map.shape[1]):
            if (contact_map[i][j] == 2) & ((j in explored) == False):
                if np.abs(j-i) in range(3,9):
                    count_alpha += 1
                    pace_helix = (np.abs(j - i) * (count_alpha - 1) + pace_helix)/count_alpha
                    explored = np.append(explored,j)
                    explored = np.append(explored,i)
            elif (contact_map[i][j] == 1) & ((j in explored) == False):
                if np.abs(j-i) > 3:
                    count_beta += 1
                    length_beta = (np.abs(j-i)*(count_beta - 1) + length_beta)/count_beta
                    explored = np.append(explored,i)
                    explored = np.append(explored,j)
    folding['pace_helix'] = pace_helix
    folding['n_carbonalpha'] = count_alpha
    percentage_alpha = count_alpha/(count_alpha+count_beta)*100
    folding['perc_alpha'] = percentage_alpha
    folding['length_beta'] = length_beta
    folding['n_carbonbeta'] = count_beta
    percentage_beta = count_beta/(count_beta+count_alpha)*100
    folding['perc_beta'] = percentage_beta
    folding['explored'] = explored
    return folding
    
    
