# This module is used to perform cumulative analysis
# on more than one trajectory object.
#
# Trajectory objects are made from module trajectory.py
#
# Here you will find functions to plot 


import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interpolate
import plotly.graph_objects as go


# WORKING
def _trajectory_antiparallel_positive(self):
    '''
    Gather data regarding antiparallel positive peptides contact from a trajectory
    
    Parameters
    -------
    self : morphoscanner.trajectory.trajectory object
    
    Returns
    -------
    x : list
        list of int in range(0, max(self.peptide_length_list)).
        where max(self.peptide_length_list) is the number of aminoacids
        in the peptide of the system with the most aminoacids.
        
    y : list
        list of float that are the timesteps of the trajectory.
        
    z : numpy.ndarray
        ndarray of shape(i,j) where i=len(y) and j=len(x).
        z(i,j) is the ratio (ap+ / tc), between [0, 1]
        where ap+ is the antiparallel positive contacts of shift x[j] at timestep y[i] 
        and tc are the total contact of that timestep.
    '''

    # Read timesteps from trajectory
    index = self.database.index
    # get timestep of each frame
    # time is in picosecond, so time/1000 for time in nanoseconds
    y = [self.universe.trajectory[i].time/1000 for i in index]
    # get the shift range
    # maximum range is the length (number of aminoacids) of the longest peptide
    x = [i for i in range(max(self.peptide_length_list))]
    # calculate total contacts per frame
    total_contact = [(self.database.iloc[i]['parallel'] + self.database.iloc[i]['antiparallel']) for i in range(len(self.database))]
    
    z = []
    # for each sampled frame
    for f_index, frame in enumerate(self.frames):

        try:
            # compute contact ratio
            f = [i/total_contact[f_index] if i!=0 else 0 for i in self.frames[frame].results.shift_profile_antiparallel_positive.values()]
        except:
            # if no contact, fill with 0
            max_shift = max(self.peptide_length_list)
            f = [0 for k in range(max_shift)]
        z.append(f)
    z = np.asarray(z)
    return x, y, z


def _trajectory_antiparallel_negative(self):
    '''
    Gather data regarding antiparallel negative peptides contact from a trajectory
    
    Parameters
    -------
    self : morphoscanner.trajectory.trajectory object
    
    Returns
    -------
    x : list
        list of int in range(0, max(self.peptide_length_list)).
        where max(self.peptide_length_list) is the number of aminoacids
        in the peptide of the system with the most aminoacids.
        
    y : list
        list of float that are the timesteps of the trajectory.
        
    z : numpy.ndarray
        ndarray of shape(i,j) where i=len(y) and j=len(x).
        z(i,j) is the ratio (ap- / tc), between [0, 1]
        where ap- is the antiparallel negative contacts of shift x[j] at timestep y[i] 
        and tc are the total contact of that timestep.
    '''

    # Read timesteps from trajectory
    index = self.database.index
    # get timestep of each frame
    # time is in picosecond, so time/1000 for time in nanoseconds
    y = [self.universe.trajectory[i].time/1000 for i in index]
    # get the shift range
    # maximum range is the length (number of aminoacids) of the longest peptide
    x = [i for i in range(max(self.peptide_length_list))]
    # calculate total contacts per frame
    total_contact = [(self.database.iloc[i]['parallel'] + self.database.iloc[i]['antiparallel']) for i in range(len(self.database))]
    
    z = []
    # for each sampled frame
    for f_index, frame in enumerate(self.frames):

        try:
            # compute contact ratio
            f = [i/total_contact[f_index] if i!=0 else 0 for i in self.frames[frame].results.shift_profile_antiparallel_negative.values()]
        except:
            # if no contact, fill with 0
            max_shift = max(self.peptide_length_list)
            f = [0 for k in range(max_shift)]
        z.append(f)
    z = np.asarray(z)
    return x, y, z

def _trajectory_parallel(self):
    '''
    Gather data regarding parallel peptides contact from a trajectory
    
    Parameters
    -------
    self : morphoscanner.trajectory.trajectory object
    
    Returns
    -------
    x : list
        list of int in range(0, max(self.peptide_length_list)).
        where max(self.peptide_length_list) is the number of aminoacids
        in the peptide of the system with the most aminoacids.
        
    y : list
        list of float that are the timesteps of the trajectory.
        
    z : numpy.ndarray
        ndarray of shape(i,j) where i=len(y) and j=len(x).
        z(i,j) is the ratio (p / tc), between [0, 1]
        where p is the parallel contacts of shift x[j] at timestep y[i] 
        and tc are the total contact of that timestep.
    '''

    # Read timesteps from trajectory
    index = self.database.index
    # get timestep of each frame
    # time is in picosecond, so time/1000 for time in nanoseconds
    y = [self.universe.trajectory[i].time/1000 for i in index]
    # get the shift range
    # maximum range is the length (number of aminoacids) of the longest peptide
    x = [i for i in range(max(self.peptide_length_list))]
    # calculate total contacts per frame
    total_contact = [(self.database.iloc[i]['parallel'] + self.database.iloc[i]['antiparallel']) for i in range(len(self.database))]

    z = []
    # for each sampled frame
    for f_index, frame in enumerate(self.frames):

        try:
            # compute contact ratio
            f = [i/total_contact[f_index] if i!=0 else 0 for i in self.frames[frame].results.shift_profile_parallel.values()]
        except:
            # if no contact, fill with 0
            max_shift = max(self.peptide_length_list)
            f = [0 for k in range(max_shift)]
        z.append(f)
    z = np.asarray(z)
    return x, y, z

def _plot_3d_antiparallel_positive_shift(x,y,z):
    fig = go.Figure(data=[go.Surface(z=z*100, x=x, y=y)])
    fig.update_layout(autosize=True,
                          scene = dict(
                        xaxis_title='AP+ Shift',
                        yaxis_title='Time (ns)',
                        zaxis_title='Contact %',
                        zaxis = dict(nticks=20, range=[0,100])),
                        title='Antiparallel Positive Shift')
    fig.show()
    return

def _plot_3d_antiparallel_negative_shift(x,y,z):
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

def _plot_3d_parallel_shift(x,y,z):
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


def plot_3d_antiparallel_positive_average(data: list):
    
    mp = [i for i in map(_trajectory_antiparallel_positive, data)]
    
    x_av = max([i[0] for i in mp])
    y_av = max([i[1] for i in mp])
    z_av = np.mean(np.asarray([i[2] for i in mp]), axis=0)
    
    _plot_3d_antiparallel_positive_shift(x_av, y_av, z_av)
    
    return
    
def plot_3d_antiparallel_negative_average(data: list):
    
    mp = [i for i in map(_trajectory_antiparallel_negative, data)]
    
    x_av = max([i[0] for i in mp])
    y_av = max([i[1] for i in mp])
    z_av = np.mean(np.asarray([i[2] for i in mp]), axis=0)
    
    _plot_3d_antiparallel_negative_shift(x_av, y_av, z_av)
    
    return
    
def plot_3d_parallel_average(data: list):
    
    mp = [i for i in map(_trajectory_parallel, data)]
    
    x_av = max([i[0] for i in mp])
    y_av = max([i[1] for i in mp])
    z_av = np.mean(np.asarray([i[2] for i in mp]), axis=0)
    
    _plot_3d_parallel_shift(x_av, y_av, z_av)
    
    return
    

def get_beta_data(self):

    index = self.database.index
    beta = [sum(i) for i in self.database['n° of peptides in macroaggregates']]
    tss_int = np.array([self.universe.trajectory[i].time/1000 for i in index]).astype(int)
    number_of_peptides = len(self.frames[0].peptides)
    x = np.linspace(tss_int.min(),tss_int.max(), tss_int.max())
    spl = interpolate.interp1d(tss_int, beta, kind ='cubic')
    beta_smooth = spl(x)
    beta_smooth_norm = (beta_smooth/number_of_peptides) * 100

    return x, beta_smooth_norm

def plot_beta_average(data:list, label: str):

    beta_data = np.array([i for i in map(get_beta_data, data)])
    beta_average = np.mean(beta_data, axis=0)[1]
    x = np.max(beta_data, axis=0)[0]

    plt.plot(x, beta_average, '-',label=label)
    plt.title('% of peptides involved in β-sheets')
    plt.ylim((0,100))
    plt.xlabel('Time (ns)')
    plt.ylabel('% of Peptides in β-sheet')
    plt.legend()
    return


def get_contacts_data(self):
    index = self.database.index
    contact = [i+e for i, e in zip(self.database['parallel'], self.database['antiparallel'])]
    antiparallel = self.database['antiparallel']
    antip_total_ratio = [anti/cont if cont != 0 else 0 for anti, cont in zip(antiparallel, contact)]
    tss_int = np.array([self.universe.trajectory[i].time/1000 for i in index]).astype(int)
    x = np.linspace(tss_int.min(),tss_int.max(), tss_int.max())
    spl = interpolate.interp1d(tss_int, antip_total_ratio, kind = 'cubic') # kind='linear'
    antip_total_ratio_smooth = spl(x)

    return x, antip_total_ratio_smooth

def plot_contacts_average(data: list, label: str):

    contacts = np.array([i for i in map(get_contacts_data, data)])
    x = np.max(contacts, axis=0)[0]
    antiparallel_contacts_average_ratio = np.mean(contacts, axis=0)[1]

    plt.plot(x, antiparallel_contacts_average_ratio,'-', label=label)
    plt.title('β-sheets alignment over time')
    plt.xlabel('Time (ns)')
    plt.ylabel('β-Sheet Organizational Index')
    plt.ylim(0,1)
    plt.legend()

    return


def get_aggregates(self):
    index = self.database.index
    aggregates = self.database['n° of macroaggreates']
    tss_int = np.array([self.universe.trajectory[i].time/1000 for i in index]).astype(int)
    x = np.linspace(tss_int.min(),tss_int.max(),tss_int.max())
    spl = interpolate.interp1d(tss_int, aggregates, kind='cubic')
    aggregates_smooth = spl(x)

    return x, aggregates_smooth

def plot_aggregates_average(data: list, label: str):
    
    aggregates_data = np.array([i for i in map(get_aggregates, data)])
    x = np.max(aggregates_data, axis=0)[0]
    y = np.mean(aggregates_data, axis=0)[1]    

    y_max = max([len(i.frames[0].peptides) for i in data])//2
    
    plt.plot(x, y,'-',label=label)
    plt.title('Aggregation Order')
    plt.xlabel('Time (ns)')
    plt.ylabel('N° of macroaggregates')
    plt.legend()
    plt.yticks([i for i in range(0,y_max+2, 2)])
    return