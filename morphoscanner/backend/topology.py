"""
Created on Thu Mar 19 20:03:20 2020

@author: lillo
"""

#from __future__ import absolute_import
import MDAnalysis as mda
from .readGro import clean_gro
import tqdm



#get a list of the number of residues of every peptide in the topology
def get_peptide_length_list(topology):
    '''
    Take the .gro file and get back a list, in which each
    element is the length in atoms of a single peptide.
    
    It relies on the number in front of the residue name
    in the .gro file.

    Parameters
    ----------
    topology : .gro file path as 'user/data/file.gro'
        The path of the .gro file on your system

    Returns
    -------
    peptide_length_list : list
        It is a list of integer, each is the number of
        residues that forms a single entry in your .gro.
        It follows the order in which you have the entry
        in the .gro file, so for example:
            
            You have:
                2 peptides of 12 residues,
                3 peptides of 8 residues
                2 seeds made from 4 peptides of 8 residues each
                
            Your output will be:
                
                [12, 12, 8, 8, 8, 32, 32]
                
        So if you have seeds in your .gro simulation data,
        each seed will be recognized as a single peptide.
                
            

    '''
    
    topology = clean_gro(topology)

    peptide_length_list = []

    temporary_list = []

    # iterate trough topology
    for residue in topology:

        # if temporary list just started, add aminoacid position in chain
        if len(temporary_list) == 0:
            temporary_list.append(int(residue[1]))

        else:
            # if position of actual residue is less than last residue
            if temporary_list[-1] > int(residue[1]):

                # append lenght of last peptide to peptide length list
                peptide_length_list.append(len(temporary_list))

                # empty temporary list
                temporary_list = []

                # append actual residue position
                temporary_list.append(int(residue[1]))

            # if position of actual residue is higher than last residue, ad current residue position
            else:
                temporary_list.append(int(residue[1]))

    # append last peptide length to length stack
    peptide_length_list.append(len(temporary_list))

    return peptide_length_list


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
        Move to memory for faster (~100x faster)
        frames coordinate retrival.
        Needs a lot of memory.

    Returns
    -------
    universe : MDAnalysis.Universe()
        DESCRIPTION.

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




# make trajectory dict

def get_coordinate_dict_from_trajectory(trj_gro, trj_xtc, peptide_length=None, start_from=0, interval=1):
    '''Parse coordinate from a .gro topology and a .xtc trajectory.
    
    Arguments:  .gro topology file path,
            
                .xtc trajectory file path,
                
                
                optional:
                    
                    peptide_length, default=None.   You can set the length of the peptide
                                                    Useful if you have to analyze simulation in which
                                                    there are premade aggregate
                                    
                    start_from, default=0.    You can chose from which frame start the counter.
                                                Useful if you are working with a simulation
                                                made of different part. Eg. If part 1 end at
                                                frame 500, you can set start_from=500 and analyze
                                                part 2 of the simulation. Use it expecially
                                                if you are sampling (interval != 1)
                                                
                    interval, default=1     Interval between sample. If you want all the frame,
                                            interval=1 (default).
                                            If you want to skip sample, this parameter let you
                                            choose the interval between 2 sample frame.
    
    
    '''

    peptides_list = get_peptide_length_list(trj_gro)

    universe = make_universe(trj_gro, trj_xtc)


    if peptide_length == None:


        n_pep = len(peptides_list)


    else:

        n_pep = sum([(e // peptide_length) for e in peptides_list])


    interval = int(interval)
    #peptide_length = peptide_length
    
    trj_dict = {}

    for index_ts, ts in tqdm.tqdm(enumerate(universe.trajectory)):

        updated_index = (index_ts + start_from)

        if (updated_index % interval) == 0:

            trj_dict[updated_index] = {}

            for peptide in range(n_pep):

                trj_dict[updated_index][peptide] = {}


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

                        for index, atom in enumerate(universe.residues[res + (counter)].atoms):


                            atom_type = str(atom).split()[2]

                            if atom_type == 'BB':

                                atom_number = (int(str(atom).split()[1].split(':')[0]) - 1)

                                coordinate = universe.atoms[atom_number].position

                                position = len(trj_dict[updated_index][peptide])


                                trj_dict[updated_index][peptide][position] = coordinate


                else:

                    for res in range(peptide_length):

                        for index, atom in enumerate(universe.residues[res + (counter)].atoms):


                            atom_type = str(atom).split()[2]


                            if atom_type == 'BB':

                                atom_number = (int(str(atom).split()[1].split(':')[0]) - 1)

                                coordinate = universe.atoms[atom_number].position

                                position = len(trj_dict[updated_index][peptide])


                                trj_dict[updated_index][peptide][position] = coordinate

                            else:
                                pass


    return trj_dict

