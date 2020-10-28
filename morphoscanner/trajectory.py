import morphoscanner
from morphoscanner import backend, trj_object
from morphoscanner.backend import distance_tensor, pattern_recognition, graph, topology

import tqdm
from timeit import default_timer as timer
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import plotly.graph_objects as go
from scipy.interpolate import interpolate


class trajectory:
    '''Class to operate on trajectory files.
    It makes an object that contain the trajectory of the simulation.
    From this object is possible to conduct analysis'''

    def __init__(self, trj_gro, trj_xtc, select = None):
        
        self.trj_gro = trj_gro
        self.trj_xtc = trj_xtc
        self.universe = topology.make_universe(self.trj_gro, self.trj_xtc)
        self.number_of_frames = len(self.universe.trajectory)
        
        if select == None:
            select = ['aminoacids']
            
        self.select = select
       
        self.peptide_length_list = topology.get_peptide_length_list(self.trj_gro, self.select)
        
        self.len_dict = topology.get_peptide_length_dict(self.peptide_length_list)
        
        print('In your trajectory there are %d frames.\n' % self.number_of_frames)

        topology.print_peptides_length(self.len_dict)
        
        return            

    
    def explore(self):
        frame = 0
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
    
    # add something to ask for threshold in main.py
    def analysis(self, frame, threshold_multiplier=1.45, device='cpu'):
        # check if threshold is given
        try:
            threshold = self.contact_threshold
        except:
            dic_0 = self.get_frame(0)
            frame_distance_0 = distance_tensor.compute_distance_and_contact_maps(dic_0, threshold=0, contacts_calculation=False, device=device)
            threshold = distance_tensor.get_median_c_alpha_distance(frame_distance_0) * threshold_multiplier
            self.contact_threshold = threshold
            print("Two nearby atoms of different peptides are contacting if the distance is lower than: %s Angstrom" % str(self.contact_threshold))
    
        #frame = frame
        print('Analyzing frame n° ', frame)
    
        frame_dict = self.get_frame(frame)
        
        start_dist = timer()
        frame_distance, frame_contact = distance_tensor.compute_distance_and_contact_maps(frame_dict, threshold=threshold, device=device)
        end_dist = timer()
        print('Time to compute distance is: ', (end_dist - start_dist))

        start_den = timer()

        frame_denoised, df = pattern_recognition.denoise_contact_maps_torch_v1(frame_contact, device=device)

        end_den = timer()
        print('Time to denoise: ', (end_den-start_den))
    
        frame_graph_full = graph.graph_v1(frame_denoised, df)
        
        subgraphs = graph.find_subgraph(frame_graph_full)  
        
        self.frames[frame].results = trj_object.trj_objects.results()
        self.frames[frame].results.distance_maps = frame_distance
        self.frames[frame].results.contact_maps = frame_contact
        self.frames[frame].results.cross_correlation = df
        self.frames[frame].results.graph = frame_graph_full
        self.frames[frame].results.subgraphs = subgraphs
        print('Finished analysis of frame n° %d' % frame)
        
        return
    
    
    def analyze_inLoop(self, threshold=None, threshold_multiplier=1.45, device='cpu'):
        
        if threshold != None:
            self.contact_threshold=threshold
        else:
            pass
        
        print('processing started...')
        start = timer()
        for frame in self.frames:
            start_an = timer()
            self.analysis(frame, threshold_multiplier=threshold_multiplier, device=device)
            end_an = timer()
            text = 'Time needed to analyze frame %d was %f seconds' % (frame, (end_an-start_an))
            print(text)

        end = timer()

        print('Total time to analyze dataset was %f seconds' % (end -start))
        return
    
    ###
    ### THESE HAVE BEEN PORTED FROM OLD TRAJECTORY TO STREAMLINE ANALYSIS OF GLICOSILATED PEPTIDES!
    ###
    
    
    def get_sense(self):

        ''' Analyze self.frames to retrieve the number of contact 
            per sense ("parallel" and "antiparallel")
        '''
        # instantiate main dict
        sense_dict = {}
        # loop trough frames
        for frame in self.frames:
            group = self.frames[frame].results.cross_correlation.groupby('sense').groups
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

                subgraph_dict[key] = morphoscanner.backend.graph.find_subgraph(self.frames[key].results.graph)

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
    
    
    def macroaggregate_sense_data(self):

        macroaggregate_sense_dict = {}

        for frame in self.frames:
            graph = self.frames[frame].results.graph
            subs = self.frames[frame].results.subgraphs
            #senses = contact_sense_in_subgraph(graph, subs)
            #sense_counter = count_sense_in_subgraph(senses)
            sense_counter = morphoscanner.backend.graph.sense_in_subgraph(graph, subs)
            macroaggregate_sense_dict[frame] = sense_counter

        self.macroaggregate_df = pd.DataFrame.from_dict(macroaggregate_sense_dict, orient='index')

        return
    
    
    def number_of_macroaggregate_per_frame(self):
        number_of_peptide = {}
        for i in self.subgraph_size_peptide:
            number_of_peptide[i] = len(self.subgraph_size_peptide[i][0])

        self.number_of_peptide_df = pd.DataFrame.from_dict(number_of_peptide, orient='index', columns=['n° of macroaggreates'])

        return
    
    
    def shift_profile(self):
        '''Calculate shift profile for the current trajectory.

        Returns
        -------
        Save data in self.frames[frame_number].results.shift....
        '''
        # instantiate main dict to store trajectory data
        shift_profile = {}
        # iterate through frames
        for frame in self.frames:
            # create a dict for each frame to store frame data
            shift_profile[frame] = {} 
            # max shift possible is the max legth of the peptides in the trajectory
            # (in number of residues)
            max_shift = max(self.peptide_length_list)
            # get cross_corr dataframe gruped by contact sense 
            a = self.frames[frame].results.cross_correlation.groupby('sense')
            # for each contact sense
            for group in a.groups:
                # create a nested dict for the contact sense
                shift_profile[frame][group] = {}
                # if sense of the contact is parallel
                # positive or negative shift is not important
                if group == 'parallel':
                    # for dataframe index (contact data)
                    for index in a.groups[group]:
                        #shift of the current contact
                        shift = int(abs(self.frames[frame].results.cross_correlation.iloc[index]['shift']))
                        #print(int(shift))
                        # if is there, add 1 to current shift counter
                        try:
                            shift_profile[frame][group][shift] += 1
                        # if is not there, just add the first contact to counter
                        except:
                            shift_profile[frame][group][shift] = 1
                        # fill with 0 shift values not found (to plot)
                        for i in range(max_shift):
                            if i not in shift_profile[frame][group]:
                                shift_profile[frame][group][i] = 0
                        # order dict by ascending keys (to plot)
                        shift_profile[frame][group] = {k[0]:k[1] for k in sorted(shift_profile[frame][group].items())}
                        # add data to frame.results
                        self.frames[frame].results.shift_profile_parallel = shift_profile[frame][group]
                # if sense of the contact is antiparallel
                # positive or negative shift is important
                # because C and N terminal of respective peptides
                # can interact differently
                if group == 'antiparallel':
                    for index in a.groups[group]:
                        shift = int(self.frames[frame].results.cross_correlation.iloc[index]['shift'])
        
                        if shift > 0:
                            shift_profile[frame][group]['negative'] = {}
                            try:
                                shift_profile[frame][group]['negative'][shift] += 1
                            except:
                                shift_profile[frame][group]['negative'][shift] = 1
                            for i in range(max_shift):
                                if i not in shift_profile[frame][group]['negative']:
                                    shift_profile[frame][group]['negative'][i] = 0
                            shift_profile[frame][group]['negative'] = {k[0]:k[1] for k in sorted(shift_profile[frame][group]['negative'].items())}
                            self.frames[frame].results.shift_profile_antiparallel_negative = shift_profile[frame][group]['negative']
        
                        if shift <= 0:
                            shift = abs(shift)
                            shift_profile[frame][group]['positive'] = {}
                            try:
                                shift_profile[frame][group]['positive'][shift] += 1
                            except:
                                shift_profile[frame][group]['positive'][shift] = 1
                            for i in range(max_shift):
                                if i not in shift_profile[frame][group]['positive']:
                                    shift_profile[frame][group]['positive'][i] = 0
                            shift_profile[frame][group]['positive'] = {k[0]:k[1] for k in sorted(shift_profile[frame][group]['positive'].items())}
                            self.frames[frame].results.shift_profile_antiparallel_positive = shift_profile[frame][group]['positive']
        return


    def get_data(self):
        self.get_sense()
        self.subgraph_length_peptide()
        self.macroaggregate_sense_data()
        self.number_of_macroaggregate_per_frame()
        self.shift_profile()
        return
    
        
    def get_database(self):
        
        self.database = pd.concat((self.subgraph_len_pep_df, self.sense_df, self.number_of_peptide_df, self.macroaggregate_df), axis=1)

        return
    
    
    ######################
    #############################
    #####################
    
    
    def plot_contacts(self, kind='cubic'):
        index = self.database.index
        contact = [i+e for i, e in zip(self.database['parallel'], self.database['antiparallel'])]
        antiparallel = self.database['antiparallel']
        antip_total_ratio = [anti/cont if cont != 0 else 0 for anti, cont in zip(antiparallel, contact)]
        tss_int = np.array([self.universe.trajectory[i].time/1000 for i in index]).astype(int)
        x = np.linspace(tss_int.min(),tss_int.max(), tss_int.max())
        spl = interpolate.interp1d(tss_int, antip_total_ratio, kind = kind)
        antip_total_ratio_smooth = spl(x)
        
        plt.plot(x, antip_total_ratio_smooth,'-')
        plt.title('β-sheets alignment over time')
        plt.xlabel('Time (ns)')
        plt.ylabel('β-Sheet Organizational Index')
    
        return
    
    
    def plot_peptides_in_beta(self, kind='cubic'):
        index = self.database.index
        beta = [sum(i) for i in self.database['n° of peptides in macroaggregates']]
        tss_int = np.array([self.universe.trajectory[i].time/1000 for i in index]).astype(int)
        number_of_peptides = len(self.frames[0].peptides)
        x = np.linspace(tss_int.min(),tss_int.max(), tss_int.max())
        spl = interpolate.interp1d(tss_int, beta, kind = kind)
        beta_smooth  = spl(x)
        beta_smooth_norm = (beta_smooth/number_of_peptides) * 100

        plt.plot(x, beta_smooth_norm, '-')
        plt.title('% of peptides involved in β-sheets')
        plt.ylim((0,100))
        plt.xlabel('Time (ns)')
        plt.ylabel('% of Peptides in β-sheet')
    
        return
    
    def plot_aggregates(self, kind='cubic'):
        index = self.database.index
        tss_int = np.array([self.universe.trajectory[i].time/1000 for i in index]).astype(int)
        aggregates = self.database['n° of macroaggreates']
        x = np.linspace(tss_int.min(),tss_int.max(), tss_int.max())
        spl = interpolate.interp1d(tss_int, aggregates, kind=kind)
        aggregates_smooth = spl(x)
        y_max = len(self.frames[0].peptides)//2

        plt.plot(x, aggregates_smooth,'-')
        plt.yticks([i for i in range(0, y_max+2, 2)])
        plt.title('Aggregation Order')
        plt.xlabel('Time (ns)')
        plt.ylabel('N° of macroaggregates')

        return
    
    def plot_shift_parallel(self, frame=None):
        try:
            f = self.frames[frame].results.shift_profile_parallel
        except:
            max_shift = max(self.peptide_length_list)
            f = {k:0 for k in range(max_shift)}
        x = [val for val in f.keys()]
        y = [k for k in f.values()]
        plt.plot(x, y)
        plt.xlabel('P Shift')
        plt.ylabel('Number of contacts')
        plt.show()
        return
    
    def plot_shift_antiparallel_positive(self, frame=None):
        try:
            f = self.frames[frame].results.shift_profile_antiparallel_positive
        except:
            max_shift = max(self.peptide_length_list)
            f = {k:0 for k in range(max_shift)}    
        x = [val for val in f.keys()]
        y = [k for k in f.values()]
        plt.plot(x, y)
        plt.xlabel('A+ Shift')
        plt.ylabel('Number of contacts')
        plt.show() 
        return
    
    def plot_shift_antiparallel_negative(self, frame=None):
        try:
            f = self.frames[frame].results.shift_profile_antiparallel_negative
        except:
            max_shift = max(self.peptide_length_list)
            f = {k:0 for k in range(max_shift)} 
            
        x = [val for val in f.keys()]
        y = [k for k in f.values()]
        plt.plot(x, y)
        plt.xlabel('A- Shift')
        plt.ylabel('Number of contacts')
        plt.show() 
        return


    def get_subgraphs_sense(self, frame):
        '''
        Retrive information about contact sense of each aggregate
        found in self.frames[frame]['subgraphs_full']
        
        Parameters
        ----------
        frame : int
            The frame of which you want to get contact sense informations.
        
        Returns
        -------
        sense_dict : dict
            A dict containing the informations about contacts, in the form:
                {'parallel' : int,
                 'antiparallel' : int,
                 'value' : str}
            The key 'value' contains the sense of the predominant contact sense,
            'parallel' or 'antiparallale',
            or the str 'equal' if both sense have the same number of contacts.
        '''
        
        # check if requested frame have been parsed
        if frame not in self.frames:
            print('Frame %d is not in the sampled frames\n' % frame)
        else:
            # check if in the frame there are aggregate
            if len(self.frames[frame].results.subgraphs) < 1:
                print('There are no aggregate in frame %d.\n' % frame)
            else:
                # if checks are passed
                # create empty dict
                sense_dict = {}
                
                # iterate subgraphs
                for index_sub, subgraph in enumerate(self.frames[frame].results.subgraphs):
                    
                    # create a new dict for each aggregate, to store contact sense information
                    sense_dict[index_sub] = {'parallel' : 0,
                                             'antiparallel' : 0,
                                             'value' : 0   }
                    # get information about contacts from database
                    #  use only peptide1 column to gather contacts one time only 
                    for index_contact, contact in enumerate(self.frames[frame].results.cross_correlation.peptide1):
                        if contact in subgraph:
                            sense = (self.frames[frame].results.cross_correlation.iloc[index_contact].sense)
                            # add 1 to the right sense counter in the sense_dict
                            sense_dict[index_sub][sense] += 1
                    # check if contacts number is equal in both senses
                    if sense_dict[index_sub]['parallel'] == sense_dict[index_sub]['antiparallel']:
                        sense_dict[index_sub]['value'] = 'equal'
                    else:
                        # if contacts are not equal, get the predominant contact sense
                        sense_dict[index_sub]['value'] = max(sense_dict[index_sub], key=sense_dict[index_sub].get)
    
                return sense_dict
    
    
    def plot_frame_aggregate(self, frame: int):
        '''
        Plot the frame with color code that identify the
        sense of the majority of contacts in an aggregate.
        
            - Grey: no contact,
            - Green: majority of parallel contacts,
            - Blue: majority of antparallel contacts,
            - Yellow: equal number of parallel and antiparallel contacts
        
        The plot can be made interactive using jupyter-notebook,
        with:
            %matplotlib notebook
        
        Parameters
        ----------
        frame : int
            The frame that you want to plot
        
        Returns
        -------
        plot
            Return a matplotlib.pyplot 3d scatter plot.
        '''
        
        # get predominant contact sense for each aggregate
        sense_dict = self.get_subgraphs_sense(frame)
        # get subgraphs
        subgraphs = self.frames[frame].results.subgraphs
        # get coordinate dict
        coordinate_dict = self.get_frame(frame)
        # make a flat (1D) list of peptide in the aggregates
        flat_subgraphs = [pep for group in subgraphs for pep in group]
        # create a color dictionary with each sense corresponding to a color
        colors = {'parallel' : 'limegreen',
                  'antiparallel' : 'b',
                  'equal' : 'y',
                  'no' : 'gray'}    
        
        # instantiate empty dict to plot aggregates
        x = {}
        y = {}
        z = {}
        # iterate through aggregates
        for index_sub, subgraph in enumerate(subgraphs):
            # create a list to gather coordinate of each aggregate's atom
            x[index_sub] = []
            y[index_sub] = []
            z[index_sub] = []
            # for each peptide in the aggregate
            for peptide in subgraph:
                # for each atom of the peptide
                for atom in coordinate_dict[peptide]:
                    # get x, y and z coordinates and save it in the correct list
                    x[index_sub].append(coordinate_dict[peptide][atom][0])
                    y[index_sub].append(coordinate_dict[peptide][atom][1])
                    z[index_sub].append(coordinate_dict[peptide][atom][2])
        
        # instantiate lists for non contacting peptides
        x_not = []
        y_not = []
        z_not = []
        # get coordinate of non contacting peptides
        for pep in coordinate_dict:
            if pep not in flat_subgraphs:
                for atom in coordinate_dict[pep]:
                    x_not.append(coordinate_dict[pep][atom][0])
                    y_not.append(coordinate_dict[pep][atom][1])
                    z_not.append(coordinate_dict[pep][atom][2])
        
        fig = plt.figure()
    
        ax = plt.axes(projection='3d')
    
        # scatter aggregates atoms
        for group in x:
    
            ax.scatter3D(x[group],y[group],z[group], color=colors[sense_dict[group]['value']])
        
        # scatter non contacting peptides atoms
        ax.scatter3D(x_not, y_not, z_not, color=colors['no'])
        
        return plt.show()
    
    
    def plot_graph(self, frame: int):
        '''
        Plot the frame graph, with visual information about
        number of contacts between peptides and sense of the contacts.
        
        Edge thickness scale with the number of contacts between two
        contacting peptides.
        
        Green edges are parallel contacts.
        Blue edges are antiparallel contacts.
        
        Parameters
        ----------
        frame : int
            The frame of which you want to plot the graph.
        
        Returns
        -------
        plot
            matplotlib.pyplot 3d scatter5 plot.
        '''
    
        graph = self.frames[frame].results.graph
        
        # Used to plot
        edges = graph.edges()
        colors = [graph[u][v][0]['color'] for u,v in edges]
        weights = [graph[u][v][0]['weight'] for u,v in edges]
        
        # output a plot
        #return nx.draw_networkx(graph, edges=edges, edge_color=colors, width=weights) #for networkx 2.4
        return nx.draw_networkx(graph, edge_color=colors, width=weights)  # for networkx 2.5


    ### Use this to plot 3d data from trajectory object
    def plot3d_parallel(self):
        # Read timestep from trajectory
        index = self.database.index
        # get timestep of each frame (in nanoseconds)
        y = [self.universe.trajectory[i].time/1000 for i in index]
        # get the shift range
        x = [i for i in range(max(self.peptide_length_list))]
        # calculate total contacts per frame
        total_contact = [(self.database.iloc[i]['parallel'] + self.database.iloc[i]['antiparallel']) for i in range(len(self.database))]
        z = []
        for f_index, frame in enumerate(self.frames):

            try:
                # compute contact ratio
                f = [i/total_contact[f_index] if i!=0 else 0 for i in self.frames[frame].results.shift_profile_parallel.values()]

            except:
                # if no contact, fill with 0
                max_shift = max(self.peptide_length_list)
                f = [0 for k in range(max_shift)]
             # append frame data
            z.append(f)
        # cast to np.array
        z = np.asarray(z, dtype=object)
        fig = go.Figure(data=[go.Surface(z=z*100, x=x, y=y)])
        fig.update_layout(autosize=True,
                          scene = dict(
                        xaxis_title='P Shift',
                        yaxis_title='Time (ns)',
                        zaxis_title='Contact %',
                        zaxis = dict(nticks=20, range=[0,100])),
                            title='Parallel Shift')
        fig.show()

        return


    def plot3d_antiparallel_negative(self):
        # Read timestep from trajectory
        index = self.database.index
        # get timestep of each frame (in nanoseconds)
        y = [self.universe.trajectory[i].time/1000 for i in index]
        # get the shift range
        x = [i for i in range(max(self.peptide_length_list))]
        # calculate total contacts per frame
        total_contact = [(self.database.iloc[i]['parallel'] + self.database.iloc[i]['antiparallel']) for i in range(len(self.database))]
        z = []
        for f_index, frame in enumerate(self.frames):

            try:
                # compute contact ratio
                f = [i/total_contact[f_index] if i!=0 else 0 for i in self.frames[frame].results.shift_profile_antiparallel_negative.values()]

            except:
                # if no contact, fill with 0
                max_shift = max(self.peptide_length_list)
                f = [0 for k in range(max_shift)]
             # append frame data
            z.append(f)
        # cast to np.array
        z = np.asarray(z, dtype=object)
        fig = go.Figure(data=[go.Surface(z=z*100, x=x, y=y)])
        fig.update_layout(autosize=True,
                          scene = dict(
                        xaxis_title='AP- Shift',
                        yaxis_title='Time (ns)',
                        zaxis_title='Contact %',
                        zaxis = dict(nticks=20, range=[0,100])),
                            title='Antiparallel Negative Shift')
        fig.show()

        return

    def plot3d_antiparallel_positive(self):
        # Read timestep from trajectory
        index = self.database.index
        # get timestep of each frame (in nanoseconds)
        y = [self.universe.trajectory[i].time/1000 for i in index]
        # get the shift range
        x = [i for i in range(max(self.peptide_length_list))]
        # calculate total contacts per frame
        total_contact = [(self.database.iloc[i]['parallel'] + self.database.iloc[i]['antiparallel']) for i in range(len(self.database))]
        z = []
        for f_index, frame in enumerate(self.frames):

            try:
                # compute contact ratio
                f = [i/total_contact[f_index] if i!=0 else 0 for i in self.frames[frame].results.shift_profile_antiparallel_positive.values()]

            except:
                # if no contact, fill with 0
                max_shift = max(self.peptide_length_list)
                f = [0 for k in range(max_shift)]
             # append frame data
            z.append(f)
        # cast to np.array
        z = np.asarray(z, dtype=object)
        fig = go.Figure(data=[go.Surface(z=z*100, x=x, y=y)])
        fig.update_layout(autosize=True,
                          scene = dict(
                        xaxis_title='AP+ Shift',
                        yaxis_title='Time (ps)',
                        zaxis_title='Contact %',
                        zaxis = dict(nticks=20, range=[0,100])),
                            title='Antiparallel Positive Shift')
        fig.show()
    
        return
    