"""
Created on Wed Nov 18 17:50:56 2020

@author: lillo
"""
import numpy as np


def contact_map_helix_torch(distance_map):
    
    t_h_sub = distance_map.clone()
    t_h_sub[t_h_sub<=4.7] = 0
    t_h_sub[(t_h_sub>4.7) & (t_h_sub<=6.5)] = 2
    t_h_sub[t_h_sub>6.5] = 3
    
    return t_h_sub


def contact_tracer_v1(contact_map):
    contact_map = contact_map.clone()
    second_type = 0
    count_second = 0 
    folding = {}
    for i in range(contact_map.shape[0]):
        for j in range(i+1, contact_map.shape[1]):
            if contact_map[i][j] == 2:
                if np.abs(j-i) == 3:
                    count_second += 1
                    second_type = (np.abs(j - i) * (count_second - 1) + second_type)/count_second

    folding['pace_alpha311'] = second_type
    folding['n_alpha311'] = count_second
    if count_second <= 2:
        return 0
    else:
        return folding
    
    