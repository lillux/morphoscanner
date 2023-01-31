"""
Created on Wed Nov 18 17:50:56 2020

@author: lillo
"""
import numpy as np
from timeit import default_timer as timer
import pandas as pd

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
    explored_alpha = []
    explored_beta = []
    for i in range(contact_map.shape[0]):
        for j in range(i+1, contact_map.shape[1]):
            if (contact_map[i][j] == 2) & ((j in explored_alpha) == False) & ((j in explored_beta) == False):
                if np.abs(j-i) in range(3,6):
                    count_alpha += 1
                    pace_helix = (np.abs(j - i) * (count_alpha - 1) + pace_helix)/count_alpha
                    explored_alpha = np.append(explored_alpha,j)
                  # if ((i in explored_alpha) == False):
                  #     explored_alpha = np.append(explored_alpha,i)
            elif (contact_map[i][j] == 1) & ((j in explored_beta) == False) & ((j in explored_alpha) == False):
                if np.abs(j-i) > 3:
                    count_beta += 1
                    length_beta = (np.abs(j-i)*(count_beta - 1) + length_beta)/count_beta
                    explored_beta = np.append(explored_beta,j)
                   # if ((i in explored_beta) == False):
                   #     explored_beta = np.append(explored_beta,i)
    folding['pace_helix'] = pace_helix
    folding['n_carbonalpha'] = count_alpha
    percentage_alpha = count_alpha/(contact_map.shape[0])*100
    folding['perc_alpha'] = percentage_alpha
    folding['length_beta'] = length_beta
    folding['n_carbonbeta'] = count_beta
    percentage_beta = count_beta/(contact_map.shape[0])*100
    folding['perc_beta'] = percentage_beta
    folding['explored_alpha'] = np.sort(explored_alpha).astype(np.int32)
    folding['explored_beta'] = np.sort(explored_beta).astype(np.int32)
    return folding


def retrieve_map(self, frame, i: int, j: int):
    ij_map = self.frames[frame].results.distance_maps[i][j]
    return ij_map


def _calculate_helix_score(self, frame):
    h_score = {}
    for peptide in range(len(self.frames[frame].peptides)):
        #print('checking peptide', peptide)
        d_map = retrieve_map(self,frame,peptide,peptide)
        #start=timer()
        h_map = contact_map_helix_torch(d_map)
        #end=timer()
        #print('contact map computed in', end-start, 'seconds')
        #start = timer()
        score = contact_tracer(h_map)
        #end = timer()
        #print('score computed in', end-start, 'seconds')
        h_score[peptide] = score
        #print('end peptide', peptide)
    return h_score

def calculate_helix_score_for_frame(self, frame: int, device='cpu'):
    h_score = _calculate_helix_score(self, frame)
    self.frames[frame].results.helix_score = pd.DataFrame(h_score).transpose()
    
    return

def helix_score(self, device='cpu'):
    '''
    Calculate alpha-helix score for each sampled timestep

    Parameters
    ----------
    device : str, optional
        The default is 'cpu'.
        Choose between 'cpu' and 'cuda'

    Returns
    -------
    None.

    '''
    for frame in self.frames:
        calculate_helix_score_for_frame(self, frame=frame, device=device)
    return


def get_max_helix_single_peptide(self, peptide:int=0):
    '''
    Find the frame in which a certain peptide has the max percentage of alpha-helix,
    and the value of the percentage, as a key:value pair.
    
    =============================================================================
     Need a morphoscanner.trajectory.trajectory() object with the analysis done,
     trajectory.helix_score() included, to work.
    =============================================================================

    Parameters
    ----------
    peptide : int, optional
        The index of the peptide to analyze.
        The default is 0.

    Returns
    -------
    dict
        a dict of which the key is the frame index,
        and the value is the percentage of alpha helix in that frame.
        
        {frame_index : percentage}

    '''
    alpha_perc = []
    for f in self.frames:
        alpha_perc.append(self.frames[f].results.helix_score.iloc[peptide]['perc_alpha'])
    max_alpha = max(alpha_perc)
    max_time = alpha_perc.index(max_alpha)
    return {max_time : max_alpha}