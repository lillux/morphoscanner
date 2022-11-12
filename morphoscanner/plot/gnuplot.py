"""
Created on Mon Oct 31 15:40:52 2022

@author: lillo

Module containing GnuPlot functionalities.
Functions works as methods of trajectory() class.
"""
import time, os, tempfile
from PyGnuplot import gp




def gnuplot_data_parallel(self):
    # Read timestep from trajectory
    index = self.database.index
    # get timestep of each frame (in nanoseconds)
    y = [self.universe.trajectory[i].time for i in index]
    # get the shift range
    x = [i for i in range(max(self.peptide_length_list))]
    # calculate total contacts per frame
    total_contact = [(self.database.iloc[i]['parallel'] + self.database.iloc[i]['antiparallel']) for i in range(len(self.database))]
    z = []
    for f_index, frame in enumerate(self.frames):
        try:
            f = [(i/total_contact[f_index])*100 if i!=0 else 0 for i in self.frames[frame].results.shift_profile_parallel.values()]
        except:
            # if no contact, fill with 0
            max_shift = max(self.peptide_length_list)
            f = [0 for k in range(max_shift)]
         # append frame data
        z.append(f)

    return x,y,z

def gnuplot_data_antiparallel_positive(self):
    # Read timestep from trajectory
    index = self.database.index
    # get timestep of each frame (in nanoseconds)
    y = [self.universe.trajectory[i].time for i in index]
    # get the shift range
    x = [i for i in range(max(self.peptide_length_list))]
    # calculate total contacts per frame
    total_contact = [(self.database.iloc[i]['parallel'] + self.database.iloc[i]['antiparallel']) for i in range(len(self.database))]
    z = []
    for f_index, frame in enumerate(self.frames):
        try:
            f = [(i/total_contact[f_index])*100 if i!=0 else 0 for i in self.frames[frame].results.shift_profile_antiparallel_positive.values()]
        except:
            # if no contact, fill with 0
            max_shift = max(self.peptide_length_list)
            f = [0 for k in range(max_shift)]
         # append frame data
        z.append(f)

    return x,y,z

def gnuplot_data_antiparallel_negative(self):
    # Read timestep from trajectory
    index = self.database.index
    # get timestep of each frame (in nanoseconds)
    y = [self.universe.trajectory[i].time for i in index]
    # get the shift range
    x = [i for i in range(max(self.peptide_length_list))]
    # calculate total contacts per frame
    total_contact = [(self.database.iloc[i]['parallel'] + self.database.iloc[i]['antiparallel']) for i in range(len(self.database))]
    z = []
    for f_index, frame in enumerate(self.frames):
        try:
            f = [(i/total_contact[f_index])*100 if i!=0 else 0 for i in self.frames[frame].results.shift_profile_antiparallel_negative.values()]
        except:
            # if no contact, fill with 0
            max_shift = max(self.peptide_length_list)
            f = [0 for k in range(max_shift)]
         # append frame data
        z.append(f)

    return x,y,z

def gnuplot_shift_surface(self, sense:str, z_max:int=100, quality:str='medium', sleep:int=2, verbose:bool=True):
    '''
    Plot a beta-sheet shift profile 3D plot in gnuplot.
    It Needs PyGnuplot installed in the system

    Parameters
    ----------
    sense : str
        is the sense of the shift that you want to plot.
        it can take one of the following values:
            'parallel'
            'antiparallel_negative'
            'antiparallel_positive'
            
    z_max : int, optional
        is the maximum value of the z axis. It start from 0,
        and finish at z_max. It is a % value,
        so the value of z_max should be 100 at maximum.
        The default is 100.
    quality : str, optional
        the quality of the saved image.
        The accepted value are:
            low',
            'medium',
            'high'
            The default is 'medium'.
    sleep : int, optional
        is the length of the pause for the interpreter, in seconds.
        This argument will be removed in future version.
        The default is 2.
    verbose : bool, optional
        if verbose is True, the function gives:
            the name of the temporary data object,
            a message confirming the finishing of the plotting
            a message confirming the deletion of the temporary file.
            The default is True.

    Raises
    ------
    ValueError
            If the value is wrong, return accepted values.
            The accepted values are ones accepted by 'sense' argument.
    Returns
    -------
    None.

    '''
    if sense == 'parallel':
        x,y,z = gnuplot_data_parallel(self)
    elif sense == 'antiparallel_negative':
        x,y,z = gnuplot_data_antiparallel_negative(self)
    elif sense == 'antiparallel_positive':
        x,y,z = gnuplot_data_antiparallel_positive(self)
    else:
        raise ValueError(f'{sense} is not a legal value. Accepted value for sense are: "paralle", "antiparallel_negative", "antiparallel_positive"')

    sense_dict={'parallel' : {'title': 'set title "Parallel Shift"',
                            'ylabel': 'set ylabel "P Shift"',
                            'output': 'set output "parallel.png"'},
               'antiparallel_negative':{'title': 'set title "Antiparallel Negative Shift"',
                            'ylabel': 'set ylabel "AP - Shift"',
                            'output': 'set output "antiparallel_neg.png"'},
               'antiparallel_positive':{'title': 'set title "Antiparallel Positive Shift"',
                            'ylabel': 'set ylabel "AP + Shift"',
                            'output': 'set output "antiparallel_pos.png"'}}
    
    quality_dict = {'low': ('640,480', '12'),
                    'medium': ('1280,960', '20'),
                    'high' : ('2560, 1920', '40')}
    # create and write on a temporary file with data to plot
    with tempfile.NamedTemporaryFile(mode='w', encoding='utf-8', delete=False) as f:
        for idy, timestep in enumerate(y):
            for idx, shift in enumerate(x):
                f.writelines(f'{timestep} {shift} {z[idy][idx]}\n')
            f.writelines('\n')
        if verbose:
            print(f.name)
        gnup = gp()
        #set title
        gnup.c(sense_dict[sense]['title'])
        # set x label (timestep)
        gnup.c(f'set xlabel "Time ({self.time_unit})"')
        # set y label (type of shift)
        gnup.c(sense_dict[sense]['ylabel'])
        # set z label (% of interaction)
        gnup.c('set zlabel "% Interaction" rotate parallel')
        # set z range
        gnup.c(f'set zrange [0:{z_max}]')
        # set pm3d mode for surface style
        gnup.c('set pm3d')
        # set opacque plot
        gnup.c('set hidden3d')
        # set plot very close to the bottom
        gnup.c('set ticslevel 0.1')
        # set view angle. TODO: make it exposed?
        gnup.c('set view 55, 230, 1, 1')
        # plot
        gnup.c(f'splot "{f.name}" w l notitle')
        # set image output
        gnup.c(f'set terminal pngcairo size {quality_dict[quality][0]} font "Verdana,{quality_dict[quality][1]}"')
        gnup.c(sense_dict[sense]['output'])
        gnup.c('replot')
    # pause interpreter, waiting for gnuplot.
    # TODO: use subprocess and wait finishing signal    
    time.sleep(sleep)
    if verbose:
        print('end')
    # remove temporary file containing the data
    os.remove(f.name)
    # exit from gnuplot
    gnup.quit()
    if verbose:
        print(f'{f.name} deleted')
    return