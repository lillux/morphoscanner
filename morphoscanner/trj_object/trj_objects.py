class frames(object):
    
    pass
    
        
class single_peptide():
    
    ''' Class that define peptides
    '''
    
    def __init__(self, seq, atom_n):
        
        self.sequence = seq
        self.atom_numbers = atom_n
        #self.frames_coordinates = frames()
        
        return
    
    
#    def get_coordinate_from_frame(self, frame, coordinates):
#        
#        name = 'frame_' + str(frame)
#
#        setattr(self.frames_coordinates, name, coordinates)
#        
#        return
    
    def get_coordinate_from_frame(self, frame, coordinates):
        
        try:
            self.frames[frame] = coordinates
        except:
            self.frames = {}
            self.frames[frame] = coordinates
        return