from timeit import default_timer as timer
from ..backend import distance_tensor, pattern_recognition, graph

        
### USING THIS ONE
# to complete
# i want to get frame and the tensor inside frame object
# then continue analysis

class frames():
    
    def __init__(self, frame):
        
        self.frame = frame
        
        pass
        

    def get_frame(self):

        a_frame = {}
        for pep in super().frames[self.frame].peptides:
            a_frame[pep] = super().frames[self.frame].peptides[pep].coordinates    
    
        pass
    
        

# Classes in dev

class single_peptide():
    
    ''' Class that define peptides
    
    '''
    
    def __init__(self, sequence, atom_n, coordinates):
        
        self.sequence = sequence
        self.atom_numbers = atom_n
        self.coordinates = coordinates
        
        return

    # do this for each peptide to gather distances
    # this is not optimized but still faster than before
#    def distance(self):
        
#        tt = distance_tensor.get_coordinate_tensor_from_dict_single(t_test.frames[0].peptides[0].coordinates)

#       self.distances = {}
#       for tens in frame_tensor:
#
#            dists[tens] = morphoscanner.backend.distance_tensor.fast_cdist(frame_tensor[tens], tt.unsqueeze(0))

#       return dists
        