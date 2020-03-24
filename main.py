# -*- coding: utf-8 -*-



import sys
import os

from morphoscanner.trajectory import trajectory


def get_gro(counter = 0):
    _gro = input('Insert the path of the .gro file: ')
    #assert os.path.exists(_gro), "I did not find the file at: "+str(_gro)
    if os.path.exists(_gro) == False:
        print("I did not find the file at %s, please retry...\n " % str(_gro))
        counter += 1
        print(counter)
        if counter >= 5:
            raise sys.exit("You are wrong boy, bye bye!") 
            #return print('You are wrong boy, bye bye!')
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
            #return print('You are wrong boy, bye bye!')
        else:
            return get_xtc(counter=counter)
    else:
        print("\n%s loaded" % _xtc)
    return _xtc



def get_value(sentence, counter = 0):
    
    value = int(input(sentence))
    #assert os.path.exists(_xtc), "I did not find the file at: "+str(_xtc)
    if type(int(value)) != int:
        print("%s is not an integer, it is of type %s...\n " % (str(value), type(value)))
        counter += 1
        print(counter)
        if counter >= 5:
            raise sys.exit("Wrong wrong boy, bye bye!") 
            #return print('You are wrong boy, bye bye!')
        else:
            return get_value(counter=counter, sentence=sentence)
    else:
        print("value is %d" % value)
    return int(value)




def get_destination_dir_and_name(counter=0):
    
    # Check for input
    destination_directory = input('Insert the path in which you want to save your output: ')
    destination_directory = os.path.normpath(destination_directory)
    if os.path.exists(destination_directory) == False:
        print("%s does not look like a proper path, please retry...\n " % str(destination_directory))
        counter += 1
        print(counter)
        if counter >= 5:
            raise sys.exit("You are wrong boy, bye bye!") 
            #return print('You are wrong boy, bye bye!')
        else:
            return get_destination_dir_and_name(counter=counter)
    else:
        print("Directory %s accepted" % destination_directory)
        
        file_name = input('Name of your output: ')
        
        extension = '.xlsx'  # TODO selectable extension (with dict?)
        
        final_path = os.path.join(destination_directory, (file_name + extension))
        
    return final_path, file_name


# script
if __name__ == '__main__':

    _gro = get_gro()
    _xtc = get_xtc()
    
    trj = trajectory(_gro, _xtc)
    
    print('Your trajectory has %d frames' % trj.number_of_frames)
    print('Yout trajectory has %d BB atoms' % trj.number_of_BB_atoms)
    
    
    #peptide_length = int(input('Set the number of aminoacids in your peptides (int): '))
    peptide_length = get_value(sentence='Set the number of aminoacids in one peptide (int): ')
    #interval = int(input('Set the interval betwen sampled frames (int): '))
    interval = get_value(sentence='Set the interval between sampled frames (int): ')
    #start_from = int(input('Set the frame from which you want to start.\n(0 if you have a single simulation. len(split1) if you are analyzing split2.\nUsed to keep track of the frame index.) (int): '))
    start_from = get_value(sentence='Set the index from which you want to start.\n(0 if you have a single simulation. 0 if you are analying split1, len(split1) if you are analyzing split2.) (int): ')
    
    output_path, file_name = get_destination_dir_and_name()
    
    
    trj.compose_database(peptide_length=peptide_length, interval=interval)
    trj.analyze_inLoop()
    trj.get_data()
    trj.get_database()
    
    
    
    trj.database.to_excel(output_path, sheet_name=file_name)

