from morphoscanner.trajectory import trajectory
from morphoscanner.backend.topology import print_peptides_length
from morphoscanner.data_acquisition.script_inputs import get_gro, get_xtc, peptide_length, get_interval, start_from, get_destination_dir_and_name


# script
if __name__ == '__main__':

    _gro = get_gro()
    _xtc = get_xtc()
    
    trj = trajectory(_gro, _xtc)
    trj.explore()
    
    interval = get_interval(sentence='\nSet the interval between sampled frames (int): ')
    trj.compose_database(sampling_interval=interval)
    
    trj.analyze_inLoop()
    trj.get_data()
    trj.get_database()
    
    output_path, file_name = get_destination_dir_and_name()
    
    trj.database.to_excel(output_path, sheet_name=file_name)

