import morphoscanner
from morphoscanner import backend
from morphoscanner.backend.check_val import isInt
import sys


class trajectory:
    '''Class to operate on trajectory files.

    It makes an object that contain the trajectory of the simulation'''

    def __init__(self, trj_gro, trj_xtc):
        
        self.trj_gro = trj_gro
        self.trj_xtc = trj_xtc
        self.universe = backend.topology.make_universe(self.trj_gro, self.trj_xtc)
        self.number_of_frames = len(self.universe.trajectory)
        self.number_of_BB_atoms = len(self.universe.select_atoms('name BB'))
       
        self.peptide_length_list = backend.topology.get_peptide_length_list(self.trj_gro)
        self.len_dict = backend.topology.get_peptide_length_dict(self.peptide_length_list)
        
        print('In your trajectory there are %d frames.\n' % self.number_of_frames)
        print('In each frame there are %d BB atoms.\n' % self.number_of_BB_atoms)
        morphoscanner.backend.topology.print_peptides_length(self.len_dict)
        
        return
        
        
    def split(self, to_split: list, split_size: list):
        '''Manually split peptide_length_list in case of seeds.
        
        Input:
            to_split: list
                list of int or ints.
                Each int refers to the length of a peptides seed
                from self.len_dict.keys() that you want to split in single peptide.
                For example if in len dict there are seeds of length 96 that you want to split,
                to_split = [96]
                
            split_size: list
                list of int or ints.
                This is the size in which you want to split your to_split seeds.
                For example if you want to split your seeds of length 96 in peptides of length 12,
                split_size = [12]
                
        Output:
            Change the original self.peptide_length_list with a new list of splitted peptides.
        
        '''
        
        splitting_dict = get_splitting_dict(to_split, split_size)
        self.peptide_length_list = get_new_peptides_length(self.peptide_length_list, splitting_dict)
        print('Splitting done.\n')
        print('"peptide_length_list" attribute has been updated with the new length.')
        
        return
    
    
    def explore(self):
        
        frame = 0
        coordinate, sequence, atom_number = get_data_from_trajectory_frame(universe=self.universe, frame=frame, peptide_length_list= self.peptide_length_list)

        self.peptide = {}
        for seq, coord, atm_n in zip(sequence, coordinate, atom_number):

            self.peptide[seq] = single_peptide(sequence.get(seq), atom_number.get(atm_n))
                
            self.peptide[seq].get_coordinate_from_frame(frame=frame, coordinates=coordinate.get(coord))
        
        print('Exploration of frame %d done.\n' % frame)
        
        return
    
    
    def compose_database(self, sampling_interval):
        
        steps = [s for s in range(self.number_of_frames) if s%sampling_interval==0 and s != 0]
        for step in tqdm.tqdm(steps):
            self.universe.trajectory[step]

            for pep in self.peptide:
                c_list = {}

                for idx, i in enumerate(self.peptide[pep].atom_numbers.values()):
                    p = self.universe.atoms[i].position
                    c_list[idx] = p

                self.peptide[pep].get_coordinate_from_frame(step, c_list)

        return
        
        
    def get_frame(self, frame):
        
        peptide_dict = {}
        for pep in self.peptide:

            peptide_dict[pep] = self.peptide[pep].frames[frame]
        
        return peptide_dict
