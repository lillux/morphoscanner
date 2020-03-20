"""
Created on Thu Mar 19 20:03:20 2020

@author: lillo
"""

#from __future__ import absolute_import
import MDAnalysis as mda
import readGro
import tqdm



#get a list of the number of residues of every peptide in the topology
def get_peptide_length_list(topology):
    
    topology = readGro.clean_gro(topology)

    peptide_lenght_list = []

    temporary_list = []

    # iterate trough topology
    for residue in topology:

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

    return peptide_lenght_list




def make_universe(trj_gro, trj_xtc):
    ''' Leverage MDAnalysis.Universe() to parse trajectory file from gromacs output.

    Intput: string: system path of gro file (topology) and
                    system path of xtc file (trajectory)
                    of the file to analyze

    return: MDAnalysis.Universe()'''

    universe = mda.Universe(trj_gro,trj_xtc)

    return universe

  


# create a dict from a Universe in which each entry is a timestep of the MD simulation
def create_trajectory_dict(universe):
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
                                                    Usefull if you hae to analyze simulation in which
                                                    there are premade aggregate
                                    
                    start_from, default=0.    You can chose from which frame start the counter.
                                                Usefull if you are working with a simulation
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

        n_pep = sum([(e//peptide_length) for e in peptides_list])


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

