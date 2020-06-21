"""
Created on Thu Mar 19 20:11:57 2020

@author: lillo
"""

import torch
#import tqdm

 # instantiate 3d tensor with shape n_peptides * n_residues * n_dimension
def get_coordinate_tensor_from_dict(coordinate_dict):
    '''Convert a frame_dict to a tensor, for parallel euclidean distance calculation.
    
    Inputs: coordinate_dict, dict.
    
    Return: torch.tensor'''

    #variables wit dict dimension
    dim0 = len(coordinate_dict)
    dim1 = len(coordinate_dict[0])
    dim2 = len(coordinate_dict[0][0])

    #initialize a 0s tensor
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
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



def get_coordinate_tensor_from_dict_multi(coordinate_dict):
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
        tensor_dict[chain] = get_coordinate_tensor_from_dict_single(coordinate_dict[chain])
    return tensor_dict



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


        
        