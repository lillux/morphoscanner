import sys
import os

from morphoscanner.trajectory import trajectory
from morphoscanner.data_acquisition import get_gro, get_xtc, peptide_length, get_interval, start_from, get_destination_dir_and_name
from morphoscanner.topology import print_peptides_length


# script
if __name__ == '__main__':

    _gro = get_gro()
    _xtc = get_xtc()
    
    trj = trajectory(_gro, _xtc)
    
    print('Your trajectory has %d frames' % trj.number_of_frames)
    print('Your trajectory has %d BB atoms\n' % trj.number_of_BB_atoms)
    print('Your simulation has the following number of peptides of a certain length:\n')
    print_peptides_length(trj.len_dict)
    
    
    peptide_length = peptide_length(sentence='Set the number of aminoacids in one peptide (int): ')
    interval = get_interval(sentence='Set the interval between sampled frames (int): ')
    start_from = start_from(sentence='Set the index from which you want to start.\n\n0 if you have a single simulation.\n0 if you are analyzing split1.\nlen(split1) if you are analyzing split2.\ninteger: ')
    
    output_path, file_name = get_destination_dir_and_name()
    
    
    trj.compose_database(peptide_length=peptide_length, interval=interval)
    trj.analyze_inLoop()
    trj.get_data()
    trj.get_database()
    
    
    
    trj.database.to_excel(output_path, sheet_name=file_name)

