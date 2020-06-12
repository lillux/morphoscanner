
"""
Created on Fri Mar 20 00:40:39 2020


    class plot():

                function to plot aggregate and peptide

@author: lillo
"""

import numpy as np
import matplotlib.pyplot as plt
import torch




# plot single peptide (with autoscaling of axes)
def plot_single_peptide(peptide_coordinate_dict, centroid=False):
    x = []
    y = []
    z = []

    for residue in peptide_coordinate_dict:
        point = peptide_coordinate_dict[residue]
        x.append(point[0])
        y.append(point[1])
        z.append(point[2])


    x = np.asarray(x)
    y = np.asarray(y)
    z = np.asarray(z)

    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.scatter3D(x,y,z, c='b')


    if centroid == True:
        median_centroid = [np.median(x), np.median(y), np.median(z)]
        ax.scatter3D(median_centroid[0], median_centroid[1], median_centroid[2], c='r')

    # Create cubic bounding box to simulate equal aspect ratio
    max_range = np.array([x.max()-x.min(), y.max()-y.min(), z.max()-z.min()]).max()
    Xb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + 0.5*(x.max()+x.min())
    Yb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten() + 0.5*(y.max()+y.min())
    Zb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten() + 0.5*(z.max()+z.min())
    # Comment or uncomment following both lines to test the fake bounding box:
    for xb, yb, zb in zip(Xb, Yb, Zb):
        ax.plot([xb], [yb], [zb], 'w')

    return plt.show()


########## PLOT PEPTIDE LIST
# plot a list of peptide point cloud in 3d space.
# The box axis have arbitrary scale dependent on the aminoacids distance
# you can select to show the centroid
def plot_peptide_list(coordinate_dict, peptide_list=None, centroid=False):
    '''
    Plot peptides from a trajectory frame.
    Using jupyter-notebook, use '%matplotlib notebook' to
    plot the points cloud in 3D.

    Parameters
    ----------
    coordinate_dict : dict
        Is the dict that contains all the coordinate
        of the atoms of a single frame.
        A single frame of the output of 
        backend.topology.get_coordinate_dict_from_trajectory 
        is a coordinate_dict.
        
    peptide_list : list, optional
        The default is None. By default all the peptides
        will be plotted.
            Is a list of int. Put here the index of the peptide
            or peptides that you want to plot.
            For example [0,2,5,24,1,6] to plot
            only these peptides.
        
    centroid : bool, optional
        The default is False.
            The centroid of a peptide can be plotted
            in red together with the selected peptide.
    Returns
    -------
    3D plot
        Return a scattered 3D plot.

    '''
    
      
    # if no peptide specified, plot all
    if peptide_list == None:
        peptide_list = [p for p in coordinate_dict]


    # if there is only a single peptide to show
    # use the single peptide function to normalize axis        
    if len(peptide_list) == 1:
        
        return plot_single_peptide(coordinate_dict[peptide_list[0]])
    
    else:
        
        x = []
        y = []
        z = []
        x_median = float
        y_median = float
        z_median = float


        for peptide in range(len(peptide_list)):
            x.append([peptide])
            y.append([peptide])
            z.append([peptide])
            for aminoacid in coordinate_dict[peptide_list[peptide]]:

                point = coordinate_dict[peptide_list[peptide]][aminoacid]
                x[peptide].append(point[0])
                y[peptide].append(point[1])
                z[peptide].append(point[2])

            del x[peptide][0]
            del y[peptide][0]
            del z[peptide][0]

        if centroid == True:

            def assemble_coordinate(axis_coordinate_list):
                median_list = []
                for coordinate_set in axis_coordinate_list:
                    median = np.median(coordinate_set)
                    median_list.append(median)
                return median_list

            x_median = assemble_coordinate(x)
            y_median = assemble_coordinate(y)
            z_median = assemble_coordinate(z)


        #%matplotlib notebook

        fig = plt.figure()

        ax = plt.axes(projection='3d')


        for pep in range(len(x)):

            ax.scatter3D(x[pep],y[pep],z[pep])

            if centroid == True:

                ax.scatter3D(x_median[pep], y_median[pep], z_median[pep], c='red')


        #return  plt.show(), [x,y,z], [x_median, y_median, z_median]         
    return plt.show()



#plot from tensor
def plot_peptide_from_tensor(coordinate_dict, peptide_list, centroid=False):

    x = []
    y = []
    z = []
    x_median = float
    y_median = float
    z_median = float


    for peptide in range(len(peptide_list)):
        x.append([peptide])
        y.append([peptide])
        z.append([peptide])
        for index, aminoacid in enumerate(coordinate_dict[peptide_list[peptide]]):

            point = coordinate_dict[peptide_list[peptide]][index]
            x[peptide].append(point[0])
            y[peptide].append(point[1])
            z[peptide].append(point[2])

        del x[peptide][0]
        del y[peptide][0]
        del z[peptide][0]

    x = torch.FloatTensor(x)
    y = torch.FloatTensor(y)
    z = torch.FloatTensor(z)

    if centroid == True:

        def assemble_coordinate(axis_coordinate_list):
            median_list = []
            for coordinate_set in axis_coordinate_list:
                median = torch.median(torch.FloatTensor(coordinate_set))
                print(median)
                median_list.append(median)
            return median_list

        x_median = assemble_coordinate(x)
        y_median = assemble_coordinate(y)
        z_median = assemble_coordinate(z)

    #%matplotlib notebook

    fig = plt.figure()

    ax = plt.axes(projection='3d')


    for pep in range(len(x)):

        # scatter points, making list from torch tensor item
        ax.scatter3D([e.item() for e in x[pep]],[e.item() for e in y[pep]],[e.item() for e in z[pep]])

        if centroid == True:

            ax.scatter3D(x_median[pep].item(), y_median[pep].item(), z_median[pep].item(), c='red')


    #return  plt.show(), [x,y,z], [x_median, y_median, z_median]         
    return plt.show()




