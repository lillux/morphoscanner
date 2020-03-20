"""
Created on Thu Mar 19 20:11:57 2020

@author: lillo
"""

import torch
import tqdm


 # instantiate 3d tensor with shape n_peptides * n_residues * n_dimension
def get_coordinate_tensor_from_dict(coordinate_dict):

    #variables wit dict dimension
    dim0 = len(coordinate_dict)
    dim1 = len(coordinate_dict[0])
    dim2 = len(coordinate_dict[0][0])

    #initialize a 0s tensor
    if torch.cuda.is_available() == True:
        
        zero = torch.zeros([dim0,dim1,dim2], dtype=torch.float32, device='cuda')

    else:
        zero = torch.zeros([dim0,dim1,dim2], dtype=torch.float32)


    for peptide in coordinate_dict:

        for aminoacid in coordinate_dict[peptide]:
            
            
            if torch.cuda.is_available() == True:
                # create torch tensor on cuda device with cordinate [x,y,z...]
                zero[peptide][aminoacid] = torch.cuda.FloatTensor(coordinate_dict[peptide][aminoacid])

            else:
                zero[peptide][aminoacid] = torch.FloatTensor(coordinate_dict[peptide][aminoacid])

                
                
    return zero




#compute euclidean norm, fast
def compute_euclidean_norm_torch(coordinate_tensor):
    '''Use matrix to compute euclidean distance dataset wise
        and return a set of distance matrix for everi couple of peptides

    ****Runs in parallel on CUDA devices.

    Argument: tensor of shape n_peptide * n_residue * number of dimension (3 for 3d)

    return: tensor of shape n_peptide * n_peptide * n_residue * n_peptide

    '''

    #create tensor of 0s with shape n_pep x n_pep * n_res + n_res
    if torch.cuda.is_available() == True:
        zero = torch.zeros((coordinate_tensor.shape[0], coordinate_tensor.shape[0], coordinate_tensor.shape[1], coordinate_tensor.shape[1]), dtype=torch.float32, device='cuda')
    else:
        zero = torch.zeros((coordinate_tensor.shape[0], coordinate_tensor.shape[0], coordinate_tensor.shape[1], coordinate_tensor.shape[1]), dtype=torch.float32)
    
    #cicle on peptide
    for index1, peptide1 in tqdm.tqdm(enumerate(coordinate_tensor)):

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
    
    if torch.cuda.is_available():
        # move to system memory and cast to numpy array
        zero = zero.cpu().numpy()

    return zero

#https://discuss.pytorch.org/t/efficient-distance-matrix-computation/9065
#https://www.dropbox.com/h?preview=Parallel+Euclidean+distance+matrix+computation+on+big+datasets.pdf      



def distance_matrix_from_2d_tensor(peptide1_tensor, peptide2_tensor):
    '''Minimal function to calculate euclidean distance between two set of points
    using quadratic expansion. Thanks to:
    https://discuss.pytorch.org/t/efficient-distance-matrix-computation/9065
    
    Input: 2d tensor of shape (n,3), 2d tensor of shape (m,3)
    
    Output: 2d tensor of shape (n,m)
            
     '''

    # calculate distance
    x_norm = torch.pow(peptide1_tensor, 2).sum(1).view(-1,1)
    y_t = torch.transpose(peptide2_tensor, 0, 1)
    y_norm = torch.pow(peptide2_tensor, 2).sum(1).view(1,-1)

    distance_map = torch.sqrt(x_norm + y_norm - 2.0 * torch.mm(peptide1_tensor, y_t))

    # put transpose of distance map in lower triangle   
    #convert nan to 0  (using this instead of torch.clamp())       
    distance_map[torch.isnan(distance_map)] = 0

    return distance_map


        
        