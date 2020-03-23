"""
Created on Thu Mar 19 20:01:35 2020

@author: lillo

"""
from .topology import make_universe, get_peptide_length_list

import pandas as pd

import torch

from .distance_tensor import distance_matrix_from_2d_tensor


def get_dataframe_from_trajectory(trj_gro, trj_xtc, peptide_length = None):
    
    '''Create pandas.DataFrame from trajectory files
    
    Arguments: str(.gro topology path),
               str(.xtc trajectory path),
               int (optional)
               
    output: pandas.DataFrame
    
    The 3rd argument (peptide_length) is optional. If you leave it empty the
    function will use the topology file to dynamically group residues into peptides,
    strictly following the topology. It adapt to different peptide lengths (aminoacids number in the peptide)
    
    
    If you do molecular dynamics simulation with peptide of fixed residues number,
    you can insert the residues number and the function will parse the residues to
    compose the peptides unsing that number.
    
    for example, you use peptides with 12 aminoacids each, but you start your simulation
    using a premade seed of 4 peptides forming a beta-sheet. The gromacs topology file will
    consider this as a 64 residues peptide. But it is actually made of 4 peptides of 12 aminoacid each.
    If you set the peptide_length to 12, it will parse the premade beta sheet 
    as 4 peptides of 12 aminoacids each.
    
    '''

    universe = make_universe(trj_gro, trj_xtc)

    #topology = trj_gro

    peptides_list = get_peptide_length_list(trj_gro)


    if peptide_length == None:


        n_pep = len(peptides_list)


    else:

        n_pep = sum([(e//peptide_length) for e in peptides_list])

    #columns_name = ['atom_number','peptide_number', 'residue_name', 'residue_position', 'coordinates']
    columns_name = ['time_step','peptide_number', 'residue_position', 'residue_name', 'atom_position', 'atom_type', 'coordinates']

    # create list for a pd.DataFrame
    # as suggested in https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.append.html
    for_pandas = []

    trj_dict = {}

    for index_ts, ts in enumerate(universe.trajectory):

        trj_dict[index_ts] = {}

        for peptide in range(n_pep):

            trj_dict[index_ts][peptide] = {}

            if peptide != 0:

                # if to check peptide_length

                if peptide_length == None:

                    counter += peptides_list[peptide - 1]

                else:
                    counter += peptide_length


            else:
                counter = 0


            if peptide_length == None:

                for res in range(peptides_list[peptide]):

                    res_name = (str(universe.residues[res + counter]).split()[1].split(',')[0])#.split(',')[0])
                    res_position = int(str(universe.residues[res + counter]).split()[2].split('>')[0])#.split(',')[0])
                    #res_id = str(res_position) + '_' + res_name

                    #print(str(universe.residues[res + counter]))

                    for index, atom in enumerate(universe.residues[res + (counter)].atoms):

                        #print(atom)

                        atom_number = (int(str(atom).split()[1].split(':')[0]) - 1)

                        atom_type = str(atom).split()[2]

                        coordinate = universe.atoms[atom_number].position

                        position = len(trj_dict[index_ts][peptide])

                        trj_dict[index_ts][peptide][position] = coordinate

                        #features = [atom_number,peptide, res_name, position, coordinate]
                        features = [index_ts, peptide, res_position, res_name, position, atom_type, coordinate]

                        for_pandas.append(features)

            else:

                for res in range(peptide_length):

                    res_name = (str(universe.residues[res + counter]).split()[1].split(',')[0])#.split(',')[0])
                    res_position = int(str(universe.residues[res + counter]).split()[2].split('>')[0]) - 1 # -1 to start index from 0
                    #res_id = str(res_position) + '_' + res_name

                    #print(str(universe.residues[res + counter]))

                    for index, atom in enumerate(universe.residues[res + (counter)].atoms):

                        #print(atom)

                        atom_number = (int(str(atom).split()[1].split(':')[0]) - 1)

                        atom_type = str(atom).split()[2]

                        coordinate = universe.atoms[atom_number].position

                        position = len(trj_dict[index_ts][peptide])

                        trj_dict[index_ts][peptide][position] = coordinate

                        #features = [atom_number,peptide, res_name, position, coordinate]
                        features = [index_ts, peptide, res_position, res_name, position, atom_type, coordinate]

                        for_pandas.append(features)



    #start = timer()
    df = pd.DataFrame(for_pandas, columns=columns_name)
    #end = timer()
    #print(end-start)
    return df



# get dataframe of bb grains for a frame ()
def get_bb(dataframe, frame):

    bb_dataframe = dataframe.groupby('time_step').get_group(frame).groupby('atom_type').get_group('BB')

    return bb_dataframe





def get_peptide_tensor_from_dataframe(dataframe, step, peptide):
    
    '''Get 2d tensor with the coordinate of every atom (or grain) of a peptide, from a dataframe
    
    
    Input:  dataframe, (the dataframe where your data are).
    
            step, (frame of the trajectory from which you want to take the coordinate)
            
            peptide, (the peptide of which you want the coordinate)
            
    Output: torch.tensor, of shape(n,3), where n is the number of atoms in the peptide
    that you are considering. And 3 are the coordinate x, y, z of the atom.
    '''

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    peptide_data = dataframe.groupby('time_step').get_group(step).groupby('peptide_number').get_group(peptide)

    peptide_tensor = torch.tensor([vector for vector in peptide_data.coordinates.values], device=device)

    return peptide_tensor






def distance_maps_from_dataframe(dataframe, time_step):
    
    '''Calculate distance maps for all the peptides in a step of the md simulation.
    
    Input:  pandas.DataFrame. Made with morphoscanner.dataframe.get_dataframe_from_trajectory()
            
            int. Timestep of the simulation
            
            
    Output: dict of tensor. len(dict)==number of peptide in the simulation
                            len(dict[i])== number of peptide in the simulation
                            tensor.shape == (number of residue in peptide1, number of residue in peptide2)
                            
    '''
    
    

    #device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    number_of_peptides = len(dataframe.groupby('time_step').get_group(time_step).groupby('peptide_number'))
    #number_of_peptides = 5

    distance_dict = {}
    # iterate trought all peptides in a frame
    for peptide1 in range(number_of_peptides):

        if peptide1 not in distance_dict.keys():

            distance_dict[peptide1] = {}

        else:
            pass
        # to assemble tensor of peptide of same size
        #number_of_atoms_in_peptide1 = len(dataframe.groupby('time_step').get_group(time_step).groupby('peptide_number').get_group(peptide1))



        peptide1_tensor = get_peptide_tensor_from_dataframe(dataframe, time_step, peptide1)

        # iterate trought peptide in the upper triangle only
        for peptide2 in range(peptide1, number_of_peptides):

            # to assemble tensor of peptide of same size
            #number_of_atoms_in_peptide2 = len(dataframe.groupby('time_step').get_group(time_step).groupby('peptide_number').get_group(peptide2))


            peptide2_tensor = get_peptide_tensor_from_dataframe(dataframe, time_step, peptide2)


            distance_map = distance_matrix_from_2d_tensor(peptide1_tensor, peptide2_tensor)
            distance_dict[peptide1][peptide2] = distance_map

            if peptide2 in distance_dict.keys():

                #distance_dict[peptide2] = {}
                distance_dict[peptide2][peptide1] = distance_map.transpose(1,0)

            else:
                distance_dict[peptide2] = {}
                distance_dict[peptide2][peptide1] = distance_map.transpose(1,0)

    return distance_dict

# SO YOU CAN CALCULATE DISTANCE MAPS
# NOW YOU HAVE TO PUT THOSE IN A DATAFRAME

#return dist

#https://discuss.pytorch.org/t/efficient-distance-matrix-computation/9065
#https://www.dropbox.com/h?preview=Parallel+Euclidean+distance+matrix+computation+on+big+datasets.pdf  


    