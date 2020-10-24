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
        

def make_universe(trj_gro, trj_xtc, in_memory=False):
    
    '''
    Parameters
    ----------
    trj_gro : string
        system path of gro file (topology).
    trj_xtc : string
        system path of xtc or trr file (trajectory).
        can be provided as a single file or
        a list of consecutive trajectory files,
        as [part1.trr, part2.trr, ...]
    in_memory : bool, optional
        The default is True.
        Move data to memory for faster (~100x faster)
        frames coordinate retrival.
        Needs more memory.

    Returns
    -------
    universe : MDAnalysis.Universe()
        
    '''

    universe = mda.Universe(trj_gro, trj_xtc, in_memory=in_memory)

    return universe


# create a dict from a Universe in which each entry is a timestep of the MD simulation
def create_trajectory_dict(universe):
    
    '''
    Parse all the universe trajectory coordinate
    and put it in a dict. It does not group peptides,
    it just get all the coordinates of BB atoms for
    each trajectory frame.

    Parameters
    ----------
    universe : mdAnalysis.Universe
        It uses mdAnalysis parsing capability to read
        the trajectory data

    Returns
    -------
    trajectory_dict : dict
        Is a dict in which each keys is an integer
        that indicate the trajectory frame, starting
        from zero.
        Values are the coordinate of all the BB atoms
        of that frame.

    '''
    
    bb = universe.select_atoms('name BB')
    trajectory_dict = {}
    
    for index, time_steps in enumerate(universe.trajectory):
        trajectory_dict[index] = bb.positions
    
    return trajectory_dict


# Do not select by using the BB nomenclature
# Use instead the aminoacids names and numbers on the first element
# and compare it with the data inside molnames

def get_data_from_trajectory_frame_v2(universe, frame, select=['aminoacids']):
    # move to frame
    universe.trajectory[frame]
    
    temporary_list = []
    pep_index = 0

    coordinate_dict = {}
    residues_dict = {}
    atom_number_dict = {}

    object_dict = {} # new

    #select = ['peptide']

    accepted_costituents = []

    for element in select:
        if element in costituents.keys():
            try:
                accepted_costituents.extend(costituents.get(element))

            except:
                accepted_costituents.append(costituents.get(element))
        else:
            raise ValueError('%s is not a valid key for morphoscanner.molnames.costituents.\n' % str(select))

    for res in universe.residues:
        if res.resname in accepted_costituents:

            res_num = res.resnum - 1 # -1 because id start from 1, but indexing start from 0

            atom = res.atoms[0] # always take the first atom of the residues (backbone)

            atom_index = atom.id - 1 # -1 because id start from 1, but indexing start from 0

            atom_coordinate = atom.position

            resname = atom.resname


            if len(temporary_list) == 0:

                temporary_list.append(res_num)

                object_dict[pep_index] = {}

                coordinate_dict[pep_index] = {}
                residues_dict[pep_index] = {}
                atom_number_dict[pep_index] = {}


            else:
                if temporary_list[-1] > res_num:

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

            coordinate_dict[pep_index][res_num] = atom_coordinate
            residues_dict[pep_index][res_num] = resname
            atom_number_dict[pep_index][res_num] = atom_index

    object_dict[pep_index] = single_peptide(residues_dict[pep_index], atom_number_dict[pep_index], coordinate_dict[pep_index])
    
    return object_dict