"""
Created on Thu Mar 19 20:03:20 2020

@author: lillo
"""

#from __future__ import absolute_import
import MDAnalysis as mda
from .readGro import clean_gro
import tqdm
from ..molnames import costituents
from .check_val import isInt
from ..trj_object.trj_objects import single_peptide


def get_peptide_length_list(path : str, select=None):
    '''
    Take the .gro file and get back a list, in which each
    element i is the number of atoms (grains) of peptide i.
    
    Only atoms specified in 'select' are selected.
    
    It relies on the number in front of the residue name
    in the .gro file and on the residue name, matching
    for name in morphoscanner.molnames.costituents['select']    

    Parameters
    ----------
    path : str
        .gro file path as 'user/data/file.gro'
        The path of the .gro file on your system.
        
   select : list of str, optional
        The default is None.
        
        select is used to select which residues to count to compose the peptide.
        It should be a list like:
            ['peptide']
        Available options are the keys retrievable with:
            morphoscanner.molnames.costituents.keys()

    Returns
    -------
    peptide_len_list : list
        A list of int. Each element is the length of that peptides,
        counting only the residues selected in 'select' argument

    '''

    if select == None:
        select = ['aminoacids']

    accepted_costituents = []

    for element in select:
        if element in costituents.keys():
            try:
                accepted_costituents.extend(costituents.get(element))

            except:
                accepted_costituents.append(costituents.get(element))

    with open(path) as gro:
        peptide_len_list = []
        temporary_list = []
        for line in gro:
            splitted = line.split()
            if len(splitted) > 1:
                # Parse residue number and name (first item of the line)
                first = (splitted[0])

                text = [i for i in first if not isInt(i)]
                text_unite = ''.join(text)

                if text_unite in accepted_costituents:

                    number = [i for i in first if isInt(i)]
                    number_unite = ''.join(number)
                    number_unite = int(number_unite)

                    if len(temporary_list) == 0:
                        temporary_list.append(number_unite)
    
                    else:
                        # get only the first atom of a residue
                        if number_unite != temporary_list[-1]:

                            if number_unite < temporary_list[-1]:
                                peptide_len_list.append(len(temporary_list))
                                temporary_list = []
                                temporary_list.append(number_unite)

                            else:
                                temporary_list.append(number_unite)
                        else:
                            pass

        peptide_len_list.append(len(temporary_list))

    return peptide_len_list


def get_peptide_length_dict(peptide_length_list):
    '''Get the number of peptides in each peptide length
    
    Input: list, the output of topology.get_peptide_length_list
    
    Output: dict, length: number of peptides of that length'''    
    
    len_dict = {}

    for i in peptide_length_list:

        if i not in len_dict:
            len_dict[i] = 1
        else:
            len_dict[i] += 1
            
    return len_dict
    
        
def print_peptides_length(len_dict):
    for key, value in len_dict.items() :
        print ('Length: %d, Peptides: %d' % (key, value))
    return
        
def get_data_from_trajectory_frame_v2(universe, frame:int, select=['aminoacids']):
    '''
    This function act as a parser for the trajectory (based on MartiniCG v2.2. molnames)

    Parameters
    ----------
    universe : MDAnalysis.Universe
        The MDAnalysis.Universe from which to parse the data
    frame : int
        The trajectory frame, or timestep, from which you want to parse the data.
    select : list(), optional
        The atom types to take in account in the data collection.
        Atom types definition are taken from `morphoscanner.molnames.costituents`.
        The default is ['aminoacids'].

    Raises
    ------
    ValueError
        ValueError if the 'select' argument is not a valid definition.

    Returns
    -------
    object_dict : dict()
        A dict of object (from morphoscanner.trj_object.trj_objects.single_peptide),
                   with the following hierarchy: dict[peptide][residue],
                   indexed from 0, in order of appearance in the system configuration file.
    '''
    # move to frame
    universe.trajectory[frame]
    # temporary_list keeps track of the residue in a peptide, to know where a peptide ends
    temporary_list = []
    pep_index = 0
    res_counter = 0
    
    coordinate_dict = {}
    residues_dict = {}
    atom_number_dict = {}

    object_dict = {} # new

    accepted_costituents = []
    
    # create a checklist for the atom types of interest
    for element in select:
        if element in costituents.keys():
            try:
                accepted_costituents.extend(costituents.get(element))

            except:
                accepted_costituents.append(costituents.get(element))
        else:
            raise ValueError('%s is not a valid key for morphoscanner.molnames.costituents.\n' % str(select))
    
    # move through the residues of the frame
    for res in universe.residues:
        if res.resname in accepted_costituents:
            # the number of the residue in the peptide
            res_num = res.resnum
            atom = res.atoms[0] # always take the first atom of the residues (backbone)
            atom_index = atom.id - 1 # -1 because id start from 1, but indexing start from 0
            atom_coordinate = atom.position
            resname = atom.resname
            # check if a new peptide is starting, and get the first C-alpha
            if len(temporary_list) == 0:
                temporary_list.append(res_num)
                
                object_dict[pep_index] = {}

                coordinate_dict[pep_index] = {}
                residues_dict[pep_index] = {}
                atom_number_dict[pep_index] = {}
            
            else:
                if temporary_list[-1] > res_num:
                    # reset the atom index counter, because a new peptide is starting
                    res_counter = 0
                    # since a new peptide is starting, save the previous peptide
                    object_dict[pep_index] = single_peptide(residues_dict[pep_index], atom_number_dict[pep_index], coordinate_dict[pep_index])
                    
                    pep_index += 1
                    
                    temporary_list = []
                    temporary_list.append(res_num)

                    object_dict[pep_index] = {}

                    coordinate_dict[pep_index] = {}
                    residues_dict[pep_index] = {}
                    atom_number_dict[pep_index] = {}
                    
                else:
                    temporary_list.append(res_num)
                    res_counter += 1

            coordinate_dict[pep_index][res_counter] = atom_coordinate
            residues_dict[pep_index][res_counter] = resname
            atom_number_dict[pep_index][res_counter] = atom_index
            
    # save the last peptide
    object_dict[pep_index] = single_peptide(residues_dict[pep_index], atom_number_dict[pep_index], coordinate_dict[pep_index])
    
    return object_dict