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
    '''
    Class to operate on trajectory files.
    
    It makes an object that contain the trajectory information of the MD simulation.
    From this object is possible to conduct the analysis.
    '''

    def __init__(self, _sys_config, _sys_traj, select = None):
        '''
        Construct the trajectory() object.

        Parameters
        ----------
        _sys_config : str
            file path of the initial configuration of the system (.gro file for GROMACS)
        _sys_traj : str
            file path of the trajectory file (or files) (.xtc or .trr in GROMACS)
        select : list(str), optional
            Choose which atoms to select and analyze. The atoms can be choosen from
            one of the morphoscanner.molnames.costituents.keys().
            
            The default is None and will select only the aminoacids alpha-carbon

        Returns
        -------
        None.

        '''
        
        # save data paths as object attribute
        self._sys_config = _sys_config
        self._sys_traj = _sys_traj
        # instantiate MDAnalysis.Universe()
        self.universe = topology.make_universe(self._sys_config, self._sys_traj)
        # save the number of frames in the trajectory as object attribute
        self.number_of_frames = len(self.universe.trajectory)
        
        # if no preference for grain selection, get aminoacids C backbone grains
        if select == None:
            select = ['aminoacids']
        # save selected grains types as object attribute
        self.select = select
        # parse the _sys_config file to get the number of selected grains for each peptide
        self.peptide_length_list = topology.get_peptide_length_list(self._sys_config, self.select)
        # save the number of aminoacids for each peptide in a dict in the form:
        #    {key:value} where:
        #       'key' is an int, a certain number of aminoacids
        #       'value' is the int number of peptides that have 'key' aminoacids
        self.len_dict = topology.get_peptide_length_dict(self.peptide_length_list)
        
        print('In your trajectory there are %d frames.\n' % self.number_of_frames)
        
        topology.print_peptides_length(self.len_dict)
        
        return           

    
    def explore(self):
        '''
        Parse the first frame of the trajectory to gather peptide sequence, coordinates 
        and index number in the MDAnalysis.Universe(). 
        '''
        # the frame to parse is the first frame of the trajectory 
        frame = 0
        # instantiate the dict() that will contain the information from all the frames in the trajectory
        self.frames = {}
        # instantiate the object 'frame'
        self.frames[frame] = trj_object.trj_objects.frames(frame)
        # instantiate 'peptide' object and fill with frame data
        self.frames[frame].peptides = backend.topology.get_data_from_trajectory_frame_v2(universe=self.universe, frame=frame, select=self.select)
        # to print when done
        print('Exploration of frame %d done.\n' % frame)

        return  
    
    
    def compose_database(self, sampling_interval=1, direct_parse=False):
        '''
        Sample the trajectory() frames to gather coordinates from the peptides.
        
        The informations retrieved about the sequence and the index number of each peptide
        are the same as the one in frame 0, already parsed with the function self.explore().

        Parameters
        ----------
        sampling_interval : int, optional
            If frame % sampling interval == 0, sample the frame.
            The default is 1.
            
        direct_parse : bool, optional
            Choose between parsing the data directly from the MDAnalysis.Universe() (True)
            or use sequence and atom_index parsede in the frame 0 to recover coordinate data from each other frames.
            The default is False.

        Returns
        -------
        None.

        '''
        # create a list that contains the trajectory step to parse
        steps = [s for s in range(self.number_of_frames) if (s % sampling_interval)==0 and (s != 0)]
        # for each step
        for step in tqdm.tqdm(steps):
            # move MDAnalysis.Universe.trajectory to the selected frame
            self.universe.trajectory[step]
            # create 'frame' object
            self.frames[step] = trj_object.trj_objects.frames(step)
            # Parse each atoms full data, slower but directly parse information from each atoms.
            if direct_parse:
                # instantiate 'peptide' object and fill with frame data.
                self.frames[step].peptides = backend.topology.get_data_from_trajectory_frame_v2(universe=self.universe, frame=step, select=self.select)
            # Parse only coordinate from each atom.
            # this is way faster, but does not work if atom index change during the MD simulation.
            else:
                self.frames[step].peptides = {}
                # for each peptide in the first frame (frame 0)
                # (number of peptides and their composition is supposed to be the same in each step)
                for pep in self.frames[0].peptides:
                    # instantiate a dict that will contain the information about ol the peptide in the frame
                    c_list = {}
                    # for each atom index parsed in frame 0
                    # (atoms are supposed to mantain their index in each step of the trajectory)
                    for idx, i in enumerate(self.frames[0].peptides[pep].atom_numbers.values()):
                        # get position of atom at atom index 'i' in frame 'step' 
                        p = self.universe.atoms[i].position
                        # add position to a dict()
                        c_list[idx] = p
                    # create peptide object, that contains the data of each grain (atom) of the peptide
                    self.frames[step].peptides[pep] = trj_object.trj_objects.single_peptide(self.frames[0].peptides[pep].sequence,self.frames[0].peptides[pep].atom_numbers,c_list)

        return
    
        
    def get_frame(self, frame):
        '''
        Get the position of all the parsed atom of a trajectory frame.

        Parameters
        ----------
        frame : int
            The trajectory frame from which you want to get the atoms coordinates from.

        Returns
        -------
        a_frame : dict
            A dict with the form:
                {peptide_index : {atom_index : [x,y,z]}}

        '''
        # instantiate the dict() that will contains the peptides coordinates
        a_frame = {}
        # for each parsed peptide
        for pep in self.frames[frame].peptides:
            # get coordinate and put them in the dict()
            a_frame[pep] = self.frames[frame].peptides[pep].coordinates

        return a_frame
    
    
    def get_peptide(self, peptide):
        '''
        Get the position of a peptide in all the parsed frames.

        Parameters
        ----------
        peptide : int
            The peptide of which you want to recover coordinates in each sampled frame of the trajectory

        Returns
        -------
        a_peptide : dict()
            A dict() with the form:
                {frame_index : {atom_index : [x,y,z]}}

        '''
        # instantiate the dict() that will contains the peptide coordinates
        a_peptide = {}
        # for each parsed frame
        for frame in self.frames:
            # get the peptide coordinates
            a_peptide[frame] = self.frames[frame].peptides[peptide].coordinates
            
        return a_peptide
    
    
    # add something to ask for threshold in main.py
    def analysis(self, frame, threshold_multiplier=1.5, device='cpu'):
        '''
        Compute analysis on a frame.

        Parameters
        ----------
        frame : int
            The frame to analyze.
            
        threshold_multiplier : float, optional
            The default is 1.5.
            threshold_multiplier is a factor used to multiply the calculated
            threshold distance for contact recognition.
            The calculate threshold distance is the median distance between all the contiguos 
            apha-carbon of each peptides aminoacid in a frame.
            
            When threshold_mmultiplier is used, the theshold is calculated as:
                
                median(median(distance (p[c], p[c+1]) for c in [0, len(p)) for each p in T[f])
                
                Where T is the set of the sampled frames
                f is a frame in T
                p is a peptide in f
                c is an alpha-carbon in p
                len(p) is the number of c in p
                distance is the euclidean distance
                median is the median value
            
            This parameter is overwritten if threshold is given.
            
        device : str, optional
            The device on which the data are saved and the analysis is computed.
            The option are:
                'cpu', to perform parallelization on cpu
                'cuda', to performe parallelization on CUDA compatible device,
                    usually Nvidia GPUs.
            
            The default is 'cpu'.

        Returns
        -------
        None.

        '''
        # check if threshold distance for contact is given
        try:
            threshold = self.contact_threshold
        # if not given, calculate the threshold
        except:
            # compute measure on sampled data
            measures = distance_tensor.sample_intrapeptide_distance(self, samples=3)
            # plot histogram
            plt.hist(measures)
            # calculate and print confidence interval and median
            low, median, high = np.percentile(measures, (2.5, 50, 97.5))
            print('95% confidence interval is between',low, 'and', high,'Angstrom.\n')
            print('median value is: %f Angstrom.\n' % median)
            # compute theshold
            threshold = median * threshold_multiplier
            print('Computed threshold is %f Angstrom.\n' % threshold)
            # save the threshold as a trajectory() object's attribute
            self.contact_threshold = threshold
            # print the threshold
            print("Two nearby atoms of different peptides are contacting if the distance is lower than: %s Angstrom.\n" % str(self.contact_threshold))
    
        # print the frame to be analyzed
        print('Analyzing frame n° ', frame)
        # get coordinate of the frame
        frame_dict = self.get_frame(frame)
        # compute pointwise distance map for each peptide
        # and compute the contact map for each peptide
        # using the threshold distance as the maximum distance counted as a contact
        start_dist = timer()
        frame_distance, frame_contact = distance_tensor.compute_distance_and_contact_maps(frame_dict, threshold=threshold, device=device)
        end_dist = timer()
        print('Time to compute distance is: ', (end_dist - start_dist), 'seconds.')
        # gather data on the contact network between peptides
        start_den = timer()
        frame_denoised, df = pattern_recognition.denoise_contact_maps_torch_v1(frame_contact, device=device)
        end_den = timer()
        print('Time to denoise: ', (end_den-start_den), 'seconds.')
        # compose a graph of the contacting peptides
        frame_graph_full = graph.graph_v1(df)
        # find isolate peptides aggregate
        subgraphs = graph.find_subgraph(frame_graph_full)  
        # save computed data in the object 'results'
        self.frames[frame].results = trj_object.trj_objects.results()
        self.frames[frame].results.distance_maps = frame_distance
        self.frames[frame].results.contact_maps = frame_contact
        self.frames[frame].results.cross_correlation = df
        self.frames[frame].results.graph = frame_graph_full
        self.frames[frame].results.subgraphs = subgraphs
        # print to confirm the end of the analysis
        print('Finished analysis of frame n° %d' % frame)
        
        return
    
    
    def analyze_inLoop(self, threshold=5, threshold_multiplier=1.5, device='cpu'):
        '''
        Compute analysis on the whole sampled dataset.
        
        Parameters
        ----------
        threshold : float, optional
            The default is 5.
            threshold is the longest distance at which two points i,j are considered in contact.
            Is unitless, the unit of measurement depends on the one used in your dataset.
            5 Angstrom is used as a default value.
            
            
            if threshold == None, the threshold is calculated as follow:
                    
                    median(median(ed (p[c], p[c+1]) for c in [0, len(p)) for each p in f) * threshold_multiplier
                    
                    Where T is the set of the sampled frames
                    f is a frame in T
                    p is a peptide in f
                    c is an alpha-carbon in p
                    len(p) is the number of c in p
                    ed is the euclidean distance
                    median is the median value
            
            
            
        threshold_multiplier : float, optional
            The default is 1.5.
            threshold_multiplier is a factor used to multiply the calculated
            threshold distance for contact recognition.
            The calculate threshold distance is the median distance between all the contiguos 
            apha-carbon median of each peptides aminoacid in a frame.
            
            When threshold_mmultiplier is used, the theshold is calculated as:
                
                median(median(distance (p[c], p[c+1]) for c in [0, len(p)) for each p in T[f])
                
                Where T is the set of the sampled frames
                f is a frame in T
                p is a peptide in f
                c is an alpha-carbon in p
                len(p) is the number of c in p
                distance is the euclidean distance
                median is the median value
            
            This parameter is not used if threshold =! None.
            
        device : str, optional
            The device on which the data are saved and the analysis is computed.
            The option are:
                'cpu', to perform parallelization on cpu
                'cuda', to performe parallelization on CUDA compatible device,
                    usually Nvidia GPUs.
            
            The default is 'cpu'.
            

        Returns
        -------
        None.
        '''
        # if threshold is given as parameter
        if threshold != None:
            # threshold for the analysis is the given threshold 
            self.contact_threshold=threshold
        else:
            pass
        # print to confirm the starting of the proces
        print('processing started...\n')
        start = timer()
        # for each sampled frame
        for frame in self.frames:
            start_an = timer()
            # compute analysis
            self.analysis(frame, threshold_multiplier=threshold_multiplier, device=device)
            end_an = timer()
            # print the time needed for the analysis
            text = 'Time needed to analyze frame %d was %f seconds.\n' % (frame, (end_an-start_an))
            print(text)

        end = timer()

        print('Total time to analyze dataset was %f seconds\n.' % (end -start))
        return
    
    ###
    ### THESE HAVE BEEN PORTED FROM OLD TRAJECTORY TO STREAMLINE ANALYSIS OF GLICOSILATED PEPTIDES!
    ###
    
    
    def get_sense(self):

        '''
        Analyze self.frames to retrieve the number of contact 
        per sense ("parallel" and "antiparallel")
                
        return: pandas.Dataframe

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
        '''
        Get information about the size of the aggregates in the trajectory
        
        Argument: self
        
        return: pandas.Dataframe
        '''
        # if there are sampled frames
        if len(self.frames) > 0:
            # instantiate dict()
            self.subgraph_size_peptide = {}
            # for each sampled frame
            for key in self.frames.keys():
                # instantiate dict()
                subgraph_dict = {}
                # find aggregates (isolated communities in the graph) and save in the dict()
                subgraph_dict[key] = morphoscanner.backend.graph.find_subgraph(self.frames[key].results.graph)
                # instantiate list()
                len_list = []
                # for each cluster isolated community)
                for i in subgraph_dict[key]:
                    # len(i) is the number of peptides in a cluster, add to the list
                    len_list.append(len(i))
                # sort cluster from bigger to smaller
                len_list.sort(reverse=True)
                # save frame data in a dict()
                self.subgraph_size_peptide[key] = [len_list]
        # save pandas.Dataframe with calculated data
        self.subgraph_len_pep_df = pd.DataFrame.from_dict(self.subgraph_size_peptide, orient='index', columns=['n° of peptides in macroaggregates'])

        return
    
    
    def macroaggregate_sense_data(self):
        '''
        Get data on the sense of the contacts between the peptides 
        '''
        # instantiate dict
        macroaggregate_sense_dict = {}
        # for each frame
        for frame in self.frames:
            # get graph
            graph = self.frames[frame].results.graph
            # get peptides clusters
            subs = self.frames[frame].results.subgraphs
            # count contact sense
            sense_counter = morphoscanner.backend.graph.sense_in_subgraph(graph, subs)
            # save data in dict()
            macroaggregate_sense_dict[frame] = sense_counter
        # save data in pandas.Dataframe
        self.macroaggregate_df = pd.DataFrame.from_dict(macroaggregate_sense_dict, orient='index')

        return
    
    
    def number_of_macroaggregate_per_frame(self):
        '''
        Calculate number of peptides clusters in each frame
        '''
        # instantiate dict()
        number_of_peptide = {}
        # for each frame data on peptides clusters size
        for i in self.subgraph_size_peptide:
            # calculate number of clusters per frame and save in a dict()
            number_of_peptide[i] = len(self.subgraph_size_peptide[i][0])
        # create pandas.Dataframe from data
        self.number_of_peptide_df = pd.DataFrame.from_dict(number_of_peptide, orient='index', columns=['n° of macroaggreates'])

        return
    
    
    def shift_profile(self):
        '''
        Calculate shift profile for the current trajectory.

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
        '''
        Single function to compute and assemble data on the whole trajectory.
        '''
        # recover data
        self.get_sense()
        self.subgraph_length_peptide()
        self.macroaggregate_sense_data()
        self.number_of_macroaggregate_per_frame()
        self.shift_profile()
        
        # put everything together
        self.database = pd.concat((self.subgraph_len_pep_df, self.sense_df, self.number_of_peptide_df, self.macroaggregate_df), axis=1)

        return
    
    
    ######################
    #############################
    #####################
    
    
    def plot_contacts(self, kind='cubic'):
        '''
        Plot the ratio `antiparallel contacts` / `total contact` for each
        sampled timestep (or frame).

        Parameters
        ----------
        kind : str, optional
            The default is 'cubic'.

        Is the kind of interpolation used to plot the data.
        The same as `scipy.interpolate.interp1()`

        Returns
        -------
        None.
        '''
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
        '''
        Plot the ratio `peptide in beta` / `total peptide` for each
        sampled timestep (or frame).

        Parameters
        ----------
        kind : str, optional
            The default is 'cubic'.

        Is the kind of interpolation used to plot the data.
        The same as `scipy.interpolate.interp1()`
        
        Returns
        -------
        None.
        '''
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
        '''
        Plot `the number of macroaggregates` or clusters for each
        sampled timestep (or frame).        

        Parameters
        ----------
        kind : str, optional
            The default is 'cubic'.

        Is the kind of interpolation used to plot the data.
        The same as `scipy.interpolate.interp1()`

        Returns
        -------
        None.
        '''
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
        colors = [graph[u][v]['color'] for u,v in edges]
        weights = [graph[u][v]['weight'] for u,v in edges]
        
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
    
    
    def plot_3d_distance_map(self, frame, i, j):
        '''
        Plot 3D distance map via plotly.
        
        Parameters
        ----------
        frame : int
            The timestep (frame index) from which you want to visualize the distance map.
        i : int
            Index of the peptide i, of which you want to visualize the map
        j : int
            Index of peptide j, of which you want to vsualize the distance map with respect to peptide i.

        Returns
        -------
        None.

        '''
        d_map = self.frames[frame].results.distance_maps[i][j]
        fig = go.Figure(data=[go.Surface(z=d_map.numpy(), x=np.linspace(0,d_map.shape[0], num=d_map.shape[0]+1), y=np.linspace(0,d_map.shape[1], num=d_map.shape[1]+1))])
        fig.update_layout(autosize=True,
                          scene = dict(
                        xaxis_title='Peptide i',
                        yaxis_title='Peptide j',
                        zaxis_title='Distance'),
                        #zaxis = dict(nticks=20, range=[0,100])),
                        title='3D distance map')
        fig.show()
        return
