# -*- coding: utf-8 -*-

'''Functions used to import gro files'''

# read .gro and output a dict with the coordinate of the file
def get_coordinate_from_gro(path):
    
    def clean_gro(path):
    
    
        # open file .gro and return a list with one element per line of the .gro file
        def read_gro(path):
            gromacs_output = open(path)
    
            gro_file = []
            for line in gromacs_output:
                gro_file.append(line)
    
    
    
            gromacs_output.close()        
    
            return gro_file



        # return string in a string with numbers
        def return_if_string(string):
            digits = []
            for i in string:
                if not i.isdigit():
                    digits.append(i)
    
            string = ''.join(digits)
    
            return string
    
    
        # return numbers in a string with numbers
        def return_if_digit(string):
            digits = []
            for i in string:
                if i.isdigit():
                    digits.append(i)
    
            string = ''.join(digits)
    
            return string


        # remove first, second and last lines from gro_file and reorder information
        # FIX OPTION TO GET ENTRY RELATED TO A LABEL (as 'bb' or 'ca')
        def clean_gro_file(gro_file):
            cleaned_gro_file = []
            for aminoacid in gro_file[2:-1]:
                splitted = aminoacid.split()
                if splitted[1] == 'BB':
                    position_in_peptide = return_if_digit(splitted[0])
                    residue = return_if_string(splitted[0])
                    index = splitted[2]
                    x = splitted[3]
                    y = splitted[4]
                    z = splitted[5]
                    cleaned_gro_file.append([index, position_in_peptide, residue, x, y, z])
            return cleaned_gro_file
    
    
        gro_file = read_gro(path)
        cleaned_gro_file = clean_gro_file(gro_file)

        return cleaned_gro_file




    # create coordinate dict from cleaned_gro_file
    def get_coordinate_dict_from_cleaned_gro(cleaned_gro_file):
        
        peptide_lenght_list = []
    
        temporary_list = []
    
        # iterate trough cleaned_gro_file
        for residue in cleaned_gro_file:
    
            # if temporary list just started, add aminoacid position in chain
            if len(temporary_list) == 0:
                temporary_list.append(int(residue[1]))
    
            else:
                # if position of actual residue is less than last residue
                if temporary_list[-1] > int(residue[1]):
    
                    # append lenght of last peptide to peptide lenght list
                    peptide_lenght_list.append(len(temporary_list))
    
                    # empty temporary list
                    temporary_list = []
    
                    # append actual residue position
                    temporary_list.append(int(residue[1]))
    
                # if position of actual residue is higher than last residue, ad current residue position
                else:
                    temporary_list.append(int(residue[1]))
    
        # append last peptide lenght to lenght stack
        peptide_lenght_list.append(len(temporary_list))
    
        # create empty dict for coordinate
        peptide_coordinate_dict = {}
    
        # create an entry in dict for every peptide in the file
        for peptide in range(len(peptide_lenght_list)):
            peptide_coordinate_dict[peptide] = {}
    
            # for every residue in lenght peptide, add coordinate x, y, z
            for residue in range(peptide_lenght_list[peptide]):
                peptide_coordinate_dict[peptide][residue] = [float(coordinate) for coordinate in cleaned_gro_file[(peptide * peptide_lenght_list[peptide])+residue][3:]]
    
        return peptide_coordinate_dict
        
    cleaned_gro = clean_gro(path)
        
    gro_coordinate_dict = get_coordinate_dict_from_cleaned_gro(cleaned_gro)
        
    return gro_coordinate_dict

