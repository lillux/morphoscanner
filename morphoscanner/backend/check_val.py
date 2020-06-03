# Test

import sys


def isInt(s):
    '''Check if s is type int and return bool.
    
    Input: object
    
    Output: bool'''
    
    try:
        return float(str(s)).is_integer()
    except:
        return False

    
    
def check_int_and_return(value):
    
    '''Check int and return value, else raise ValueError and print object type
    
    Input = object
    
    Output = int'''

    if isInt(value):

        return int(value)

    else:
        raise ValueError("%s is not an integer, it is of type %s...\n " % (str(value), type(value))) 
