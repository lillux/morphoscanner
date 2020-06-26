from timeit import default_timer as timer
from ..backend import distance_tensor, pattern_recognition, graph

        
class single_peptide():
    
    ''' Class that define peptides
    '''
    
    def __init__(self, seq, atom_n):
        
        self.sequence = seq
        self.atom_numbers = atom_n
        
        return
    

    def get_coordinate_from_frame(self, frame, coordinates):
        
        try:
            self.frames[frame] = coordinates
        except:
            self.frames = {}
            self.frames[frame] = coordinates
        return
    

    def get_data(self, data_name, data):
        
        setattr(self, data_name, data)
        
    
    

class results():
    
    def __init__(self, frame):
        self.frame = frame
        pass

    def get_data(self, data_name, data):
        
        setattr(self, data_name, data)