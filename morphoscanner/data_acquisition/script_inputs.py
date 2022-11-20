import os
import sys
from ..backend.check_val import isInt

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


def check_input_multiple_int_recursive_with_sentence(sentence, limit=5):
    '''
    Ask for input and evaluate if input are int. Can take multiple input.

    Parameters
    ----------
    sentence : str
        The sentence that will be printed on screen while asking for the input.

    limit : int, optional
        The default is 5.
        Number of times to ask for input in case of wrong input

    Raises
    ------
    sys
        Exit if input is wrong 'limit' times

    Returns
    -------
    list
        Return a list of int.

    '''
    
    value = input(sentence)
    input_list = value.split()
    
    if len(input_list) == 0:
        limit -= 1
        
        if limit == 0:
            raise sys.exit("Too many empty inputs. Closing...")
        else:
            print('%d trial left.\n'
                  'You forgot to insert a value...please retry.' % limit)
            return check_input_multiple_int_recursive_with_sentence(sentence=sentence, limit=limit)
    
    else:
        va_list = []
        for val in input_list:
            
            if isInt(val):
                va_list.append(int(val))
                
            else:
                limit -= 1
                print('%d trial left.' % (limit))
                if limit == 0:
                    raise sys.exit("%s is not an integer, it is of type %s...\nClosing... " % (str(val), type(val))) 
                else:
                    print("%s is not an integer, it is of type %s...\n " % (str(val), type(val)))
                    return check_input_multiple_int_recursive_with_sentence(sentence=sentence, limit=limit)
        
        return va_list
    

def peptide_length(sentence, counter = 0):
    
    value = input(sentence)
    
    if value == '':
        value = None
        print('The .gro topology file is set as reference for the analysis')
        return value
    
    else:
        if isInt(value):
            
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


def ask_for_splitting(limit=5):
    '''
    Ask if the user wants to split the peptides.
    
    Used as part of the requests made to split peptides seeds in .gro file

    Parameters
    ----------
    limit : int, optional
        The default is 5.
        Number of times to ask for input in case of wrong input


    Raises
    ------
    sys
        Exit if input is wrong 'limit' times

    Returns
    -------
    bool
        Return True if the users wants to split, else False.

    '''

    answer = input("Do you want to split peptides? Write 'yes' or 'no': ")
    
    if answer not in {'n','no','y','yes'}:
        print('This is not a valid answer, please write yes or no.\n'
            '%d trial left.' % limit)
        limit -= 1
        if limit == 0:
            raise sys.exit('Too many wrong inputs. Closing...')
        else:
            return ask_for_splitting(limit=limit)

    elif answer in {'n', 'no'}:
        print('The .gro topology file is set as reference for the analysis')
        return False

    elif answer in {'y', 'yes'}:
        return True

    
def check_for_compatibility(list1, list2):
    if len(list1) == len(list2):
        
        for e1, e2 in zip(list1, list2):
            if e1 % e2 != 0:
                print('%d is not multiple of %d' % (e1,e2))
                return False
        else:
            return True
    
    else:
        raise ValueError('Your lists are of different len! list1 len = %d, list2 len = %d.' % (len(list1), len(list2)))


def get_splitting_dict(to_split, split_size):
    
    splitting_dict = {}
    
    for length, split_dim in zip(to_split, split_size):
        #check for divisibility
        if (length % split_dim) == 0:
            splitting_dict[length] = split_dim
    
    return splitting_dict

def get_new_peptides_length(peptide_length_list, splitting_dict):
    new_peptide_list = []
    for pep_length in peptide_length_list:
        if pep_length in splitting_dict.keys():
            new_size = splitting_dict[pep_length]
            new_peptide_list.extend([splitting_dict[pep_length] for p in range((pep_length//new_size))])
        else:
            new_peptide_list.append(pep_length)
    
    return new_peptide_list