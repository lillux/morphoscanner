"""
Created on Thu Mar 19 20:11:57 2020

@author: lillo
"""

import torch



 # instantiate 3d tensor with shape n_peptides * n_residues * n_dimension
def get_coordinate_tensor_from_dict(coordinate_dict, device='cpu'):
    '''
    Convert a frame_dict to a torch.tensor,
    for parallel euclidean distance calculation.
    
    Works only if all the peptides in the coordinate_dict 
    have the same number of aminoacids.

    Parameters
    ----------
    coordinate_dict : dict
        dict in the form {peptide_index : {atom_index : [x,y,z]}}
        
    device : str, optional
        Choose between 'cpu' or 'gpu'. The default is 'cpu'.
        It is the device in which the tensor is instantiated.
        
    Returns
    -------
    torch.tensor
        DESCRIPTION.

    '''

    #variables with dict dimensions
    dim0 = len(coordinate_dict)
    dim1 = len(coordinate_dict[0])
    dim2 = len(coordinate_dict[0][0])

    #initialize a 0s tensor
    #device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    zero = torch.zeros([dim0,dim1,dim2], dtype=torch.float32, device=device)

    for peptide in coordinate_dict:

        for aminoacid in coordinate_dict[peptide]:
            
            zero[peptide][aminoacid] = torch.tensor(coordinate_dict[peptide][aminoacid], device=device)
                
    return zero

def get_coordinate_tensor_from_dict_single(coordinate_dict, device='cpu'):
    '''
        Convert a coordinate_dict to a torch.tensor, for parallel euclidean distance calculation.
        Works on dict in the form {atom_key : [x, y, z]}

    Parameters
    ----------
    coordinate_dict : dict
        Is the coordinate_dict in the form {key : [x, y, z]}.
        It also works for N-dimensional points.

    Returns
    -------
    zero : torch.tensor
        Returns a torch.tensor of shape n x m
        'n'  are the keys in coordinate_dict al len(coordinate_dict)
        'm' is the number of dimensions of your data points
        
        It save on gpu if torch.cuda.is_available(), else on cpu
        If you want to move your data on cpu, e.g. for visualization,
        you need to output_tensor.cpu()
    '''
    
    #variables with dict dimension
    dim0 = len(coordinate_dict)
    first_key = [k for k in coordinate_dict.keys()][0]
    dim1 = len(coordinate_dict[first_key])

    #initialize a 0s tensor
    #device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    zero = torch.zeros([dim0,dim1], dtype=torch.float32, device=device)

    for index, peptide in enumerate(coordinate_dict):
            
        zero[index] = torch.tensor(coordinate_dict[peptide], device=device)
                
    return zero



def get_coordinate_tensor_from_dict_multi(coordinate_dict, device='cpu'):
    '''
    Generate tensor from multichain coordinate dict.
    Your coordinate_dict is in the form:
        
        {chain : {atom : [x, y, z] }}

    Parameters
    ----------
    coordinate_dict : dict
        Your coordinate_dict.
        It is in the form:
        {chain : {atom : [x, y, z] }}.

    Returns
    -------
    tensor_dict : dict
        It is a dict of tensor, one tensor per chain.

    '''
    tensor_dict = {}
    for chain in coordinate_dict:
        tensor_dict[chain] = get_coordinate_tensor_from_dict_single(coordinate_dict[chain], device=device)
    return tensor_dict


# group tensors of same size in a single tensor, and put them in a dict
def cat_tensor_for_size(tensor_dict):
    '''Group 2D tensors of same shape in a 3D tensor,
        and then put in a dict.    

    Parameters
    ----------
    tensor_dict : dict
        output of get_coordinate_tensor_from_dict_multi

    Returns
    -------
    container_tensor : dict
        a dict for a key for each tensor size, of which value is a 3D tensor
        that contains all the 2D tensor of same shape.
        
    peptide_index_in_tensor : dict
        a dict to track peptide index, in the form:
            {new_index: original_index}

    '''
    
    container_tensor = {}
    peptide_index_in_tensor = {}
    
    for i in tensor_dict:
        actual_tensor_len = len(tensor_dict[i])


        if actual_tensor_len in container_tensor.keys():
            container_tensor[actual_tensor_len] = torch.cat((container_tensor[actual_tensor_len],tensor_dict[i].unsqueeze(0)))
            peptide_index_in_tensor[actual_tensor_len][len(peptide_index_in_tensor[actual_tensor_len])] = i
        else:
            container_tensor[actual_tensor_len] = tensor_dict[i].unsqueeze(0)
            peptide_index_in_tensor[actual_tensor_len] = {}
            peptide_index_in_tensor[actual_tensor_len][len(peptide_index_in_tensor[actual_tensor_len])] = i


    return container_tensor, peptide_index_in_tensor


def compute_distance_between_each_peptide(coordinate_dict):
    '''Parallel euclidean distance matrix calculation for frame coordinate dict.
    Compute in parallel for frame containing peptide of different length.

    Parameters
    ----------
    coordinate_dict : dict
        Coordinate of a frame, in the form {peptide_index: {atom_index: atom_coordinate}}.

    Returns
    -------
    distance_maps_dict : dict
        A dict containing the computed distance maps, in the form:
            {peptide_i : {peptide_j : distance_maps_ij}}.

    '''
    
    coordinate_tensor = get_coordinate_tensor_from_dict_multi(coordinate_dict)
    group_tensor, order_tensor = cat_tensor_for_size(coordinate_tensor)
    
    distance_maps_dict = {}
    for tensor in coordinate_tensor:
        distance_maps_dict[tensor] = {}
        for tensor_group in group_tensor:
            index_dict = order_tensor[tensor_group]
            distance = fast_cdist(coordinate_tensor[tensor],group_tensor[tensor_group])
            for map_index, m in enumerate(distance):
                real_index = index_dict[map_index]
                distance_maps_dict[tensor][real_index] = m

    return distance_maps_dict


#compute euclidean norm, fast
def compute_euclidean_norm_torch(coordinate_tensor, device='cpu'):
    '''Use matrix to compute euclidean distance dataset wise
        and return a set of distance matrix for everi couple of peptides

    ****Runs in parallel on CUDA devices.

    Argument: tensor of shape n_peptide * n_residue * number of dimension (3 for 3d)

    return: tensor of shape n_peptide * n_peptide * n_residue * n_residue

    '''

    #create tensor of 0s with shape n_pep x n_pep * n_res + n_res
    #device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    zero = torch.zeros((coordinate_tensor.shape[0], coordinate_tensor.shape[0], coordinate_tensor.shape[1], coordinate_tensor.shape[1]), dtype=torch.float32, device = device)
    
    #cicle on peptide
    for index1, peptide1 in enumerate(coordinate_tensor):

        #cicle on peptide (upper triangle + diagonal)
        for index2 in range(index1, coordinate_tensor.shape[0]):
        #for index2 in range(coordinate_tens.shape[0]):

            #coordinate vector
            peptide2 = coordinate_tensor[index2]

            x_norm = torch.pow(peptide1, 2).sum(1).view(-1,1)
            y_t = torch.transpose(peptide2, 0, 1)
            y_norm = torch.pow(peptide2, 2).sum(1).view(1,-1)

            dist = torch.sqrt(x_norm + y_norm - 2.0 * torch.mm(peptide1, y_t))

            #dist = x_norm + y_norm - 2.0 * torch.mm(peptide1, y_t)
            #fine = torch.clamp(dist, 0.0, np.inf) #should be there, but is not working somehow

            # add distance map in the right position of the 0s tensor
            zero[index1][index2] = dist

            # if mesuring between different peptides
            if index1 != index2:
                # put transpose of distance map in lower triangle
                zero[index2][index1] = dist.transpose(1,0)

    #convert nan to 0  (using this instead of torch.clamp())       
    zero[torch.isnan(zero)] = 0
    
    #if device=='cuda' and torch.cuda.is_available():
    if device == 'cuda':
        # move to system memory and cast to numpy array
        zero = zero.cpu().numpy()
        
    else:
        zero = zero.numpy()

    return zero

#https://discuss.pytorch.org/t/efficient-distance-matrix-computation/9065
#https://hal.archives-ouvertes.fr/hal-02047514     



def distance_matrix_from_2d_tensor(peptide1_tensor, peptide2_tensor=None, device='cpu'):
    '''
    Minimal function to calculate euclidean distance between two set of points
    using quadratic expansion. Thanks to:
            https://discuss.pytorch.org/t/efficient-distance-matrix-computation/9065
            https://github.com/pytorch/pytorch/pull/25799
            https://github.com/pytorch/pytorch/issues/15253
    
    Parameters
    ----------
    peptide1_tensor : torch.tensor
        torch.tensor of shape n x d.
        
    peptide2_tensor : torch.tensor, optional
        The default is None.
        torch.tensor for which you want to calculate te distance from peptide1_tensor
        If None, this will be a copy of peptide_tensor1
        for pairwise distance matrix calculation
        shape m x p
        
    device : str, optional
        The default is 'cpu'.
        Is the device on which to compute the calculation
        You can set it to 'cuda' if you have an Nvidia GPU and CUDA driver installed

    Returns
    -------
    distance_map : torch.tensor
        shape n x p
        tensor with the distances data
        If you use 'cuda' device, you need to move output to main memory, with:

        distance_map.cpu()

    '''

    if peptide2_tensor == None:
        peptide2_tensor = peptide1_tensor

    # calculate distance
    x_norm = torch.pow(peptide1_tensor, 2).sum(1).view(-1,1)
    y_t = torch.transpose(peptide2_tensor, 0, 1)
    y_norm = torch.pow(peptide2_tensor, 2).sum(1).view(1,-1)
    
    distance_map = torch.sqrt(x_norm + y_norm - 2.0 * torch.mm(peptide1_tensor, y_t))
    
    # convert nan to 0  (using this instead of torch.clamp())       
    distance_map[torch.isnan(distance_map)] = 0
    
    # if you are calculating pointwise distance a single tensor
    # main diagonal is 0, to fix stability errors
    if peptide1_tensor is peptide2_tensor:
        distance_map = distance_map.fill_diagonal_(0)
    
    return distance_map


## This works if also to multiply a tensor with a matrix via broadcasting
## Taken from a discussion on the web
## TO PUT REFERENCE AND CITATION
## From Gaussian Process
def fast_cdist(x1, x2):
    '''Euclidean distance calculation between tensors.
    
    Parameters
    ----------
    x1 : torch.tensor
        Your first tensor to multiply.
        dims: a x b x c
        
    x2 : torch.tensor
        Your second tensor to multiply.
        dims: a x n x c
        You have to broadcast if it is a matrix
        all dims should be equal between tensors, except -2

    Returns
    -------
    res : torch.tensor
        the euclidean distance between the two tensors
        dims: a x b x n

    '''
    adjustment = x1.mean(-2, keepdim=True)
    x1 = x1 - adjustment
    x2 = x2 - adjustment  # x1 and x2 should be identical in all dims except -2 at this point

    # Compute squared distance matrix using quadratic expansion
    x1_norm = x1.pow(2).sum(dim=-1, keepdim=True)
    x1_pad = torch.ones_like(x1_norm)
    x2_norm = x2.pow(2).sum(dim=-1, keepdim=True)
    x2_pad = torch.ones_like(x2_norm)
    x1_ = torch.cat([-2. * x1, x1_norm, x1_pad], dim=-1)
    x2_ = torch.cat([x2, x2_pad, x2_norm], dim=-1)
    res = x1_.matmul(x2_.transpose(-2, -1))

    # Compute square root to have actual distance
    # Zero out negative values
    res = res.sqrt()
    res[torch.isnan(res)] = 0
    return res


### WORKING
def compute_distance_and_contact_maps(coordinate_dict, threshold=None, contacts_calculation=True, device='cpu'):

    # instantiate tensor
    coordinate_tensor = get_coordinate_tensor_from_dict_multi(coordinate_dict, device=device)
    # group tensor for size
    group_tensor, order_tensor = cat_tensor_for_size(coordinate_tensor)

    distance_maps_dict = {}
    contact_maps_dict = {}

    # move through tensor
    for tensor in coordinate_tensor:
        distance_maps_dict[tensor] = {}
        contact_maps_dict[tensor] = {}
        # move through tensor size
        for tensor_group in group_tensor:
            index_dict = order_tensor[tensor_group]

            #calculate distance
            distance = fast_cdist(coordinate_tensor[tensor],group_tensor[tensor_group])

            #save data maps in dictionary
            for map_index, m in enumerate(distance):
                real_index = index_dict[map_index]
                distance_maps_dict[tensor][real_index] = m
            
            if contacts_calculation == True:
                
                #crate mask for contact threshold
                contact_mask = distance.new_full(distance.shape, fill_value=threshold)
                #calculate contacts
                contact = distance - contact_mask
                #clean masks
                contact[contact > 0] = 0
                contact[contact < 0] = 1
                
                #save contact maps in dict 
                for cont_index, c in enumerate(contact):
                    real_index = index_dict[cont_index]
                    contact_maps_dict[tensor][real_index] = c
    
    if contacts_calculation == True:
        return distance_maps_dict, contact_maps_dict
    else:
        return distance_maps_dict


def get_median_c_alpha_distance(distance_maps):
    '''Calculate the median distance between the consecutive
    C-alpha of all peptides of a frame's distance maps.
    This function is used to compute a contact distance threshold that
    will be used to compute the contact_maps.

    Parameters
    ----------
    distance_maps : dict
        A dict containing the distance maps of a trajectory frame's peptides.

    Returns
    -------
    threshold : torch.float32
        Is the median distance between all of the consecutive C-alpha

    '''
    # instantiate an empty list
    median_list = []
    # move through peptide (key of the dict)
    for row in distance_maps:
        # get the distance map of the peptide with itself
        self_distance_map = distance_maps[row][row]
        # get the diagonal +1 (that is the diagonal containing the
        # distance between atom[i] and atom[i+1], that are consecutive
        # c-alpha.
        intrapep_distance_median = torch.median(torch.diag(self_distance_map,1))
        # calculate the median of that peptide's c-alpha distance
        median_list.append(intrapep_distance_median)
    # calculate the median distance for all the peptides
    threshold = torch.median(torch.tensor(median_list))
    return threshold
