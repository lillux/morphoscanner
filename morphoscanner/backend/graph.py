"""
    class graph():

        
        Function to assemble and analyze graph, using NetworkX
        
        https://networkx.github.io/

@author: lillo
"""


from .math_utilities import get_col_vector, get_row_vector

import numpy as np   
import networkx as nx    
    

# graph clustering
def nx_graph_search(denoised_dict, minimum_contacts = 3):
    
    ''' Graph clustering of peptides in the aggregate.
    
    Input: denoised contact maps dict
    
    return: networkx.MultiGraph
    
    #######
    
    Search for group of minimum 3 peptides (beta_sheet),
    joined at least with 'minimum_contacts'. default = 3
    
    '''
    # Instantiate graph
    graph = nx.MultiGraph()
    
    #iter on peptide
    for peptide_1 in denoised_dict:
        
        #iter on peptide joined to peptide1
        for peptide_2 in denoised_dict[peptide_1]:
            #retrieve contact_map of peptide1 and peptide2
            array_1 = denoised_dict[peptide_1][peptide_2]
            
            #iter on peptide joined to peptide2
            for peptide_3 in denoised_dict[peptide_2]:
                if peptide_3 != peptide_1:
                    #retrieve contact_map of peptide1 and peptide2
                    array_2 = denoised_dict[peptide_2][peptide_3]
                    
                    # create row and column vector from contact maps
                    vect_1 = get_row_vector(array_1)
                    vect_2 = get_col_vector(array_2)
                    
                    #check for contact
                    contacts = np.dot(vect_1, vect_2)
                    
                    #add edge only if there are enough contacts
                    if contacts >= minimum_contacts:

                        graph.add_edge(peptide_1, peptide_2)

                        graph.add_edge(peptide_2, peptide_3)

    return graph

#A novel graph clustering algorithm based on discrete-time quantum random walk
#S.G. Roya, A. Chakrabarti

# working with networkX
 
# when you add_edge, nodes are created if they are not there
# you can put info in edge (as distance, n of contacts, contact map)


# create a full graph
def full_graph(denoised_dict):
    ''' Create a full graph of all the peptides in the frame.

    Every peptide is a node in the graph.
    Edges join contacting peptides.
    Edges have attribute 'length' that gives you the number of contact between the peptides

    Useful for peptides behaviour analysis during molecular dynamics

    Arguments: denoised contact maps dict
    return: networkx.MultiGraph

    '''
    
    # Instantiate graph
    graph = nx.MultiGraph()

    for peptide_1 in denoised_dict:
        for peptide_2 in denoised_dict[peptide_1]:

            array_1 = denoised_dict[peptide_1][peptide_2]

            graph.add_node(peptide_1)
            graph.add_node(peptide_2)

            number_of_contacts = array_1.sum()

            if number_of_contacts >= 1:

                graph.add_edge(peptide_1, peptide_2, length = number_of_contacts)

    return graph


## THIS IS WORKING in object trajectory

def graph_v1(denoised_dict, df):
    ''' Create a full graph of all the peptides in the frame.

    Every peptide is a node in the graph.
    Edges join contacting peptides.
    Edges have attribute:   'length', that gives the number of contact between the peptides
                            'sense', that gives the sense of the contact

    Useful to analyze peptides aggregation behaviour analysis during molecular dynamics

    Arguments: denoised contact maps dict
    return: networkx.MultiGraph

    '''
    colors = {'parallel' : 'g',
              'antiparallel' : 'b'}
    
    # Instantiate graph
    graph = nx.MultiGraph()

    # iter dataframe (second output of denoise_contact_maps)
    for group in df.index:
        
        # peptides index assignment
        peptide_1 = df.iloc[group]['peptide1']
        peptide_2 = df.iloc[group]['peptide2']

        # take array from denoised contact map using index peptides as index
        array_1 = denoised_dict[peptide_1][peptide_2]
        
        # add peptides nodes in the graph 
        graph.add_node(peptide_1)
        graph.add_node(peptide_2)

        #check number of contacts between the peptides
        number_of_contacts = array_1.sum()

        # if there is a contact between the peptides
        # ad an edge and data attributes to edge
        if number_of_contacts >= 1:

            sense = df.iloc[group]['sense']

            graph.add_edge(peptide_1, peptide_2, weight=number_of_contacts, color=colors[sense])

    return graph


#FIND SUBGRAPH
def find_subgraph_in_order(graph):
    '''
    Find subgraph of joined peptides that have no node in common,
    starting from a node (peptide) with degree==1 (first or last node of the subgraph).
    The node are in order from start to end of the beta-sheet.
    
    DOES NOT WORK IF EACH PEPTIDE IS TOUCHING MORE THAN ONE OTHER PEPTIDE
    
    so if adjacency > 1, or degree > 1 for each node in the graph, it does not work.
    
    Argument: NetworkX MultiGraph

    Return: list of subgraph ordered from one end to the other

    '''

    subgraph_list = []

    for node in graph:

        # don't explore node that are already in subgraph_list
        if node not in set(nod for nod_list in subgraph_list for nod in nod_list):

            # tree is the list of nodes joined to node, starting from node
            # using depht first search
            tree = [e for e in nx.algorithms.traversal.depth_first_search.dfs_tree(graph, node)]

            # check if the first node of the tree has adjiacency == 1
            # so it checks if it is the first or last node of the subgraph
            if len(graph[tree[0]]) == 1:

                if len(subgraph_list) == 0:
                    subgraph_list.append(tree)

                else:
                    # use generator to check if the tree is already in the subgraph
                    if set(tree) not in (set(i) for i in subgraph_list):
                        subgraph_list.append(tree)

    return subgraph_list


def find_subgraph(graph):
    '''
    Find subgraph of joined peptides that have no node in common with other subgraph.
    It start to search from always from peptide 0, and from that does depth first search.
    The peptide of the subgraph that are touching are in consequential order in the sublist

    Argument: NetworkX MultiGraph

    Return: list of subgraph

    '''

    subgraph_list = []

    for node in graph:

        # don't explore node that are already in subgraph_list
        if node not in set(nod for nod_list in subgraph_list for nod in nod_list):

            # tree is the list of nodes joined to node, starting from node
            # using depht first search
            tree = [e for e in nx.algorithms.traversal.depth_first_search.dfs_tree(graph, node)]

            # check if the first node of the tree has adjiacency == 1
            # so it checks if it is the first or last node of the subgraph
            # DOES NOT FIND ANYTHING IF THERE ARE ALWAYS MORE THAN 1 CONTACT BETWEEN PEPTIDES
            #if len(graph[tree[0]]) == 1:

            if len(subgraph_list) == 0:
                subgraph_list.append(tree)

            else:
                # use generator to check if the tree is already in the subgraph
                if set(tree) not in (set(i) for i in subgraph_list):
                    subgraph_list.append(tree)
    
    subgraph_list.sort(key=len, reverse=True)

    return subgraph_list


def get_not_in_subgraph(coordinate_dict, subgraph):
    '''Get peptide alone or in cluster of 2 in a frame.
    
    Input: coordinate_dict, output of 'find_subgraph' function            
    
    Output: list
    
    #####
    
    Basically this function gives you all the node left out
    from the 'find_subgraph' function.
    The output depends on the graph creation function:
    
    nx_graph_search: get all the node with 0 and 1 neighbour
                    if you leave the 'minimum_contact' parameter
                    of the nx_graph_search function to the default value of 3.
                    
    graph_v1 :      get node without neighbours
    
    '''
    # one line function # don't use clever shits my friend
    #out = [e for e in coordinate_dict if e not in [a for i in subgraph for a in i]]

    out = []
    
    # get a list with all the node in subgraph
    subgraph = [a for i in subgraph for a in i]

    # iter on all element and get the one that are not in subgraph
    for e in coordinate_dict:

        if e not in subgraph:

            out.append(e)

    return out


def subgraph_length(aggregate):
    '''Get information about the size of the aggregates in the trajectory
    
    Work in object, fix to omogenizr with other graph in/out

    Argument: aggregate

    return: dict, keys = frame number,
                  value = a sorted list (big to small) of the aggregate size in that frame.
    '''

    subgraph_len_dict = {}
    subgraph_dict = {}
    
    for key in aggregate.frames.keys():

        subgraph_dict[key] = find_subgraph(aggregate.frames[key]['frame_graph_full'])

        len_list = []

        for i in subgraph_dict[key]:

            len_list.append(len(i))

        len_list.sort(reverse=True)

        subgraph_len_dict[key] = len_list

        return subgraph_len_dict
    
    
def contact_sense_in_subgraph(graph, subgraph_list):
    '''Get information about the orientation of the contact between contacting node in the graph.
    
    "subgraph_list" is similar to an "n_bunch" in NetworkX definition,
    but works if the node in the n_bunch are joined by an edge.'''

    subgraph_sense_dict = {}

    for index, subgraph in enumerate(subgraph_list):

        subgraph_sense_dict[index] = {}

        sense_list = []

        for i in graph.edges(subgraph):

            #print(i[0], i[1])
            sense = graph[i[0]][i[1]][0]['sense']
            sense_list.append(sense)

            subgraph_sense_dict[index] = sense_list

    return subgraph_sense_dict
    
      
def sense_in_subgraph(graph, subgraphs_list):
    '''Get number of contact per sense in subgraph.
    
    Arguments: NetworkX.MultiGraph, as output from morphoscanner.graph.graph_v1.
    
    Output = dict, with number of parallels and antiparallels contact per aggregate'''                   

    subgraphs = contact_sense_in_subgraph(graph, subgraphs_list)

    sense_counter_dict = {}

    for index, subgraph in enumerate(subgraphs):

        sense_counter_dict[index] = {}

        parallel = 0
        antiparallel = 0

        for sense in subgraphs[index]:

            if sense == 'parallel':
                parallel += 1

            elif sense == 'antiparallel':
                antiparallel += 1

        sense_counter_dict[index] = { 'parallel' : parallel,
                                       'antiparallel' : antiparallel}

    return sense_counter_dict


