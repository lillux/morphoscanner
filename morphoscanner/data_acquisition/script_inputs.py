import os
import sys

def get_gro(counter = 0):
    _gro = input('Insert the path of the .gro file: ')
    #assert os.path.exists(_gro), "I did not find the file at: "+str(_gro)
    if os.path.exists(_gro) == False:
        print("I did not find the file at %s, please retry...\n " % str(_gro))
        counter += 1
        print(counter)
        if counter >= 5:
            raise sys.exit("You are wrong boy, bye bye!") 
        else:
            return get_gro(counter=counter)
    else:
        print("\n%s loaded" % _gro)
    return _gro


def get_xtc(counter = 0):
    _xtc = input('Insert the path of the .xtc file: ')
    #assert os.path.exists(_xtc), "I did not find the file at: "+str(_xtc)
    if os.path.exists(_xtc) == False:
        print("I did not find the file at %s, please retry...\n " % str(_xtc))
        counter += 1
        print(counter)
        if counter >= 5:
            raise sys.exit("Wrong again boy, bye bye!") 
        else:
            return get_xtc(counter=counter)
    else:
        print("\n%s loaded" % _xtc)
    return _xtc



def peptide_length(sentence, counter = 0):
    
    value = input(sentence)
    
    if value == '':
        value = None
        print('The .gro topology file is set as reference for the analysis')
        return value
    
    else:
        if str.isdigit(value):
            
            return int(value)
        
        else:
            counter += 1
            print(counter)
            if counter >= 5:
                raise sys.exit("%s is not an integer, it is of type %s...\n " % (str(value), type(value))) 
            else:
                print("%s is not an integer, it is of type %s...\n " % (str(value), type(value)))
                return peptide_length(counter=counter, sentence=sentence)
            
            

def get_interval(sentence, counter = 0):

    value = input(sentence)
    
    if value == '':
        value = 1
        print('The interval was set to 1, each frame of the trajectory will be sampled')
        return value
    
    else:
        if str.isdigit(value):
            
            return int(value)
        
        else:
            counter += 1
            print(counter)
            if counter >= 5:
                raise sys.exit("%s is not an integer, it is of type %s...\n " % (str(value), type(value))) 
            else:
                print("%s is not an integer, it is of type %s...\n " % (str(value), type(value)))
                return get_interval(counter=counter, sentence=sentence)


def start_from(sentence, counter = 0):
    
    value = input(sentence)
    
    if value == '':
        value = 0
        print('Your analysis starts from frame index: %d' % value)
        return value
    
    else:
        if str.isdigit(value):
            value = int(value)
            print('Your analysis starts from frame index: %d' % value)
            return value
        
        else:
            counter += 1
            print(counter)
            if counter >= 5:
                raise sys.exit("%s is not an integer, it is of type %s...\n " % (str(value), type(value))) 
            else:
                print("%s is not an integer, it is of type %s...\n " % (str(value), type(value)))
                return start_from(counter=counter, sentence=sentence)



def get_destination_dir_and_name(counter=0):
    
    # Check for input
    destination_directory = input('Insert the path in which you want to save your output: ')
    destination_directory = os.path.normpath(destination_directory)
    if os.path.isdir(destination_directory) == False:
        print("%s does not look like an existing directory, please retry...\n " % str(destination_directory))
        counter += 1
        print(counter)
        if counter >= 5:
            raise sys.exit("You are wrong boy, bye bye!") 
        else:
            return get_destination_dir_and_name(counter=counter)
    else:
        print("Directory %s accepted" % destination_directory)
        
        file_name = input('Name of your output: ')
        
        extension = '.xlsx'  # TODO selectable extension (with dict?)
        
        final_path = os.path.join(destination_directory, (file_name + extension))
        
    return final_path, file_name
