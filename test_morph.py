# -*- coding: utf-8 -*-

import sys
import os

from morphoscanner.trajectory import trajectory

# part 1
#prod_xtc = '/home/lillo/TesiCNTE/from_cluster/prod/prod_part1/prod.xtc'   # laptop
#prod_gro = '/home/lillo/TesiCNTE/from_cluster/prod/prod_part1/min.gro'    # laptop

# part 2
#prod1_xtc = '/home/lillo/TesiCNTE/from_cluster/prod/prod_part2/prod-compl.xtc' #laptop
#prod1_gro = '/home/lillo/TesiCNTE/from_cluster/prod/prod_part2/prod-compl.gro' #laptop


_gro = input('Insert the path of the .gro file: ')
assert os.path.exists(_gro), "I did not find the file at: "+str(_gro)
print("\n%s loaded" % _gro)

_xtc = input('Insert the path of the .xtc file: ')
assert os.path.exists(_xtc), "I did not find the file at: "+str(_xtc)
print("\n%s loaded" % _xtc)

trj = trajectory(_gro, _xtc)



#prod.compose_database(peptide_length=18, interval=2000)
#prod.analyze_inLoop()
#prod.get_data()
#prod.get_database()

#prod.database.to_excel('/home/lillo/Documents/Code/export_test/t1.xls', sheet_name='t1')
