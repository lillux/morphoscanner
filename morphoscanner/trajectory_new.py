import morphoscanner
from morphoscanner import backend, data_acquisition, trj_object
from morphoscanner.backend import distance_tensor, pattern_recognition, graph, topology
from morphoscanner.backend.check_val import isInt

import tqdm
from timeit import default_timer as timer
import sys



class trajectory:
    '''Class to operate on trajectory files.
    It makes an object that contain the trajectory of the simulation.
    From this object is possible to conduct analysis'''

    def __init__(self, trj_gro, trj_xtc, select = None):
        
        self.trj_gro = trj_gro
        self.trj_xtc = trj_xtc
        self.universe = backend.topology.make_universe(self.trj_gro, self.trj_xtc)
        self.number_of_frames = len(self.universe.trajectory)
        
        if select == None:
            select = ['aminoacids']
            
        self.select = select
       
        self.peptide_length_list = backend.topology.get_peptide_length_list(self.trj_gro, self.select)
        
        self.len_dict = backend.topology.get_peptide_length_dict(self.peptide_length_list)
        
        print('In your trajectory there are %d frames.\n' % self.number_of_frames)

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
        
        splitting_dict = data_acquisition.script_inputs.get_splitting_dict(to_split, split_size)
        self.peptide_length_list = data_acquisition.script_inputs.get_new_peptides_length(self.peptide_length_list, splitting_dict)
        print('Splitting done.\n')
        print('"peptide_length_list" attribute has been updated with the new length.')
        
        return    
    
    def explore(self, frame=0): # you can change the frame number if you want to manually explore other frames
    
        self.frames = {}
        self.frames[frame] = trj_object.trj_objects.frames(frame)
        self.frames[frame].peptides = backend.topology.get_data_from_trajectory_frame_v2(universe=self.universe, frame=frame, select=self.select)
        print('Exploration of frame %d done.\n' % frame)

        return    
    
    def compose_database(self, sampling_interval=1):
        
        steps = [s for s in range(self.number_of_frames) if (s % sampling_interval)==0 and (s != 0)]
        for step in tqdm.tqdm(steps):
            self.universe.trajectory[step]
            self.frames[step] = trj_object.trj_objects.frames(step)
            self.frames[step].peptides = {}
            for pep in self.frames[0].peptides:
                c_list = {}

                for idx, i in enumerate(self.frames[0].peptides[pep].atom_numbers.values()):
                    p = self.universe.atoms[i].position
                    c_list[idx] = p

                self.frames[step].peptides[pep] = trj_object.trj_objects.single_peptide(self.frames[0].peptides[pep].sequence,self.frames[0].peptides[pep].atom_numbers,c_list)

        return
        
    def get_frame(self, frame):
        
        a_frame = {}

        for pep in self.frames[frame].peptides:
            a_frame[pep] = self.frames[frame].peptides[pep].coordinates

        return a_frame
    
    def get_peptide(self, peptide):
    
        a_peptide = {}
        for frame in self.frames:
            
            a_peptide[frame] = self.frames[frame].peptides[peptide].coordinates
            
        return a_peptide
    
    def analysis(self, frame):
    
        #frame = frame
        print('Analyzing frame n° ', frame)
    
        frame_dict = self.get_frame(frame)
    
        frame_tensor = distance_tensor.get_coordinate_tensor_from_dict(frame_dict)
    
        start_dist = timer()
        frame_distance_maps = distance_tensor.compute_euclidean_norm_torch(frame_tensor)
        end_dist = timer()
        print('Time to compute distance is: ', (end_dist - start_dist))
    
        start_contc = timer()
        frame_contact = pattern_recognition.compute_contact_maps_as_array(frame_distance_maps)
        end_contc = timer()
        print('Time to compute contact is: ', (end_contc - start_contc))
    
        start_den = timer()
        frame_denoised, df = pattern_recognition.denoise_contact_maps(frame_contact)
        end_den = timer()
        print('Time to denoise: ', (end_den-start_den))
    
        frame_graph_full = graph.graph_v1(frame_denoised, df)
        
        subgraphs = graph.find_subgraph(frame_graph_full)  
        
        self.frames[frame].results = trj_object.trj_objects.results()
        self.frames[frame].results.cross_correlation = df
        self.frames[frame].results.graph = frame_graph_full
        self.frames[frame].results.subgraphs = subgraphs
        print('Finished analysis of frame n° %d' % frame)
        
        return
    
        
    def analyze_inLoop(self):
        
        print('processing started...')
        start = timer()
        for frame in self.frames:
            start_an = timer()
            self.analysis(frame)
            end_an = timer()
            text = 'Time needed to analyze frame %d was %f seconds' % (frame, (end_an-start_an))
            print(text)

        end = timer()


        print('Total time to analyze dataset was %f seconds' % (end -start))
        return
    