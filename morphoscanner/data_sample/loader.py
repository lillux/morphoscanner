"""
Created on Sun Nov 20 01:38:24 2022

@author: lillo
"""
import os
import pathlib
import importlib_resources
import inspect
import importlib
import sys



class dataloader:
    '''
    Load data for testing
    '''
    
    def __init__(self, lib):
       # save lib path as a str
        self.current_path = os.path.dirname(sys.argv[0])
       # save .gro path as an os agnostic path object, then cast to str
        self.gro_simple_str = str(pathlib.Path(self.current_path).joinpath('simple_homogeneous', 'system.gro'))
       # save .xtc path as an os agnostic path object, then cast to str
        self.xtc_simple_str = str(pathlib.Path(self.current_path).joinpath('simple_homogeneous', 'trajectory.xtc'))
        return
