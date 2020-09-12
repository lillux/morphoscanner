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
        

#    def get_frame(self):
#
#        a_frame = {}
#        for pep in super().frames[self.frame].peptides:
#            a_frame[pep] = super().frames[self.frame].peptides[pep].coordinates    
#    
#        pass
    

# save results in this
class results():
    
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

