"""
Created on Thu Mar 19 12:38:50 2020

@author: lillo
"""


from . import backend
from backend.topology import get_peptide_length_dict
from timeit import default_timer as timer
import pandas as pd

class trajectory:

    '''Class to operate on trajectory files.

    It makes an object that contain the trajectory of the simulation'''



    def __init__(self, trj_gro, trj_xtc):
        
        

        self.trj_gro = trj_gro
        self.trj_xtc = trj_xtc
        self.universe = backend.topology.make_universe(self.trj_gro, self.trj_xtc)
        self.number_of_frames = len(self.universe.trajectory)
        self.number_of_BB_atoms = len(self.universe.select_atoms('name BB'))
        self.frames = {}
        
        peptide_length_list = backend.topology.get_peptide_length_list(self.trj_gro)
        self.len_dict = get_peptide_length_dict(peptide_length_list)



    def compose_database(self, peptide_length=None, start_from=0, interval=1):

        self.peptide_length = peptide_length
        self.start_from = start_from
        self.interval = interval

        self.data = backend.topology.get_coordinate_dict_from_trajectory(self.trj_gro, self.trj_xtc, peptide_length=self.peptide_length, start_from=self.start_from, interval=self.interval)
        self.sampled_frames = [key for key in self.data.keys()]





    def analysis(self, frame):

        # WHY len(frame_denoised) is len(frame_dict)-1 ???????

        self.frame = frame
        print('Analyzing frame n° ', self.frame)

        self.frame_dict = self.data[self.frame]

        self.frame_tensor = backend.distance_tensor.get_coordinate_tensor_from_dict(self.frame_dict)

        start_dist = timer()
        self.frame_distance_maps = backend.distance_tensor.compute_euclidean_norm_torch(self.frame_tensor)
        end_dist = timer()
        print('Time to compute distance is: ', (end_dist - start_dist))
        
        start_contc = timer()
        self.frame_contact = backend.pattern_recognition.compute_contact_maps_as_array(self.frame_distance_maps)
        end_contc = timer()
        print('Time to compute contact is: ', (end_contc - start_contc))
        
        start_den = timer()
        self.frame_denoised, self.df = backend.pattern_recognition.denoise_contact_maps(self.frame_contact)
        end_den = timer()
        print('Time to denoise: ', (end_den-start_den))
        
        self.frame_graph = backend.graph.nx_graph_search(self.frame_denoised)

        self.frame_graph_full = backend.graph.graph_v1(self.frame_denoised, self.df)

        self.subgraphs = backend.graph.find_subgraph(self.frame_graph_full)

        if self.frame not in self.frames:

            self.frames[self.frame] = {'frame_dict': self.frame_dict,
                                       'frame_denoised': self.frame_denoised,
                                          'frame_data' : self.df,
                                          'frame_graph' : self.frame_graph,
                                          'frame_graph_full' : self.frame_graph_full,
                                          'subgraphs_full' : self.subgraphs}            


    def analyze_inLoop(self):

        if hasattr(self, 'sampled_frames'):
            print('processing started...')
            start = timer()
            for frame in self.sampled_frames:
                start_an = timer()
                self.analysis(frame)
                end_an = timer()
                text= 'Time needed to analyze frame %d is %f seconds' % (frame, (end_an-start_an))
                print(text)

            end = timer()


            print('Total time to analyze dataset is %f seconds' % (end -start))

        else:
            print('You have to compose the database before analyze it.')

            print('Use "compose_database" attribute to make a database first!')

        return

    def get_sense(self):
        
        ''' Analyze self.frames to retrieve the number of contact 
            per sense ("parallel" and "antiparallel")
                
        '''
        
        # instantiate main dict
        sense_dict = {}
        
        # loop trough frames
        for frame in self.frames:
            
            group = self.frames[frame]['frame_data'].groupby('sense').groups
            
            # check for antiparallel key in the frame_data
            if 'antiparallel' in group:
            
                # get number of antiparallel contacts
                antiparallel = len(group['antiparallel'])
            
            else:
                antiparallel = 0
                
                
            # check for parallel key in the frame_data
            if 'parallel' in group:
                        
                # get number of parallel contacts
                parallel = len(group['parallel'])
        
            else:
                parallel = 0
            
            # add frame data to main dict
            sense_dict[frame] = {  'parallel' : parallel,
                               'antiparallel' : antiparallel}

        # at the end convert dict to pandas.DataFrame
        self.sense_df = pd.DataFrame.from_dict(sense_dict, orient='index')


        return



    # TODO add support to retrieve peptides index of peptides in macroaggregate

    def subgraph_length_peptide(self):
        '''Get information about the size of the aggregates in the trajectory

        Argument: aggregate

        return: dict, keys = frame number,
                      value = a sorted list (big to small) of the aggregate size in that frame


        '''

        if len(self.frames) > 0:


            self.subgraph_size_peptide = {}

            for key in self.frames.keys():

                subgraph_dict = {}

                subgraph_dict[key] = backend.graph.find_subgraph(self.frames[key]['frame_graph_full'])

                len_list = []

                for i in subgraph_dict[key]:

                    len_list.append(len(i))

                len_list.sort(reverse=True)

                self.subgraph_size_peptide[key] = [len_list]


        self.subgraph_len_pep_df = pd.DataFrame.from_dict(self.subgraph_size_peptide, orient='index', columns=['n° of peptides in macroaggregates'])

        #else:
         #   print('You have to analyze one or more frame before analyze the results.')
         #   print('Use "Analyze" or "AnalyzeInLoop" on the dataset first!')

        return



    def subgraph_length_seed(self):
        '''Get information about the size of the aggregates in the trajectory

        Argument: aggregate

        return: dict, keys = frame number,
                      value = a sorted list (big to small) of the aggregate size in that frame


        '''
        if len(self.frames) > 0:

            self.subgraph_size_seed = {}

            for key in self.frames.keys():

                subgraph_dict = {}

                subgraph_dict[key] = backend.graph.find_subgraph(self.frames[key]['frame_graph'])

                len_list = []

                for i in subgraph_dict[key]:

                    len_list.append(len(i))

                len_list.sort(reverse=True)

                self.subgraph_size_seed[key] = [len_list]

            #return self.seed_subgraph_len

        #else:
           # print('You have to analyze one or more frame before analyze the results.')
            #print('Use "Analyze" or "AnalyzeInLoop" on the dataset first!'')

        return




    def macroaggregate_sense_data(self):

        macroaggregate_sense_dict = {}

        for frame in self.frames:
            graph = self.frames[frame]['frame_graph_full']
            subs = self.frames[frame]['subgraphs_full']
            #senses = contact_sense_in_subgraph(graph, subs)
            #sense_counter = count_sense_in_subgraph(senses)
            sense_counter = backend.graph.sense_in_subgraph(graph, subs)
            macroaggregate_sense_dict[frame] = sense_counter

        self.macroaggregate_df = pd.DataFrame.from_dict(macroaggregate_sense_dict, orient='index')

        return
    
    
    
    def number_of_macroaggregate_per_frame(self):
        number_of_peptide = {}
        for i in self.subgraph_size_peptide:
            number_of_peptide[i] = len(self.subgraph_size_peptide[i][0])

        self.number_of_peptide_df = pd.DataFrame.from_dict(number_of_peptide, orient='index', columns=['n° of macroaggreates'])

        return

    
    
    


    def get_data(self):
        self.get_sense()
        self.subgraph_length_peptide()
        self.macroaggregate_sense_data()
        self.number_of_macroaggregate_per_frame()
        
        return
    
    
    
    
    def get_database(self):
        
        self.database = pd.concat((self.subgraph_len_pep_df, self.sense_df, self.number_of_peptide_df, self.macroaggregate_df), axis=1)

        return
        
























