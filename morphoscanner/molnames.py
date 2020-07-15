# Use this function to get name from Martini itp atoms names

def get_molnames(path):
    name_list = []
    flag = False
    with open(path) as martini_amino:
        for line in martini_amino:
            if flag:
                name_list.append(line.split()[0])
                flag = False

            if len(line.split()) > 1:
                actual_line = line.split()

                actual_line = [i.split(';') for i in actual_line]
                if 'molname' in actual_line[1] or 'molname' in actual_line[0]:
                    flag = True
    return name_list



#from martini_v2.2_aminoacids.itp
aminoacids_list =   ['GLY',
                     'ALA',
                     'CYS',
                     'VAL',
                     'LEU',
                     'ILE',
                     'MET',
                     'PRO',
                     'ASN',
                     'GLN',
                     'ASP',
                     'ASP0',
                     'GLU',
                     'GLU0',
                     'THR',
                     'SER',
                     'LYS',
                     'LYS0',
                     'ARG',
                     'ARG0',
                     'HIS',
                     'HISH',
                     'PHE',
                     'TYR',
                     'TRP',
                     'BTNr4']

#from martini_v2.2.itp
water_list = ['W', 'WF']

#from martini_v2.0_solvents.itp
solvents_list =     ['BUT',
                    'POP',
                    'OCT',
                    'DEC',
                    'HD',
                    'OD',
                    'ODE',
                    'ODT',
                    'BENZ',
                    'CHEX',
                    'BOL',
                    'POL',
                    'EOL',
                    'OCO',
                    'PRX',
                    'CLF',
                    'MES',
                    'ETH',
                    'PON',
                    'PAM',
                    'ACH',
                    'ACE',
                    'TO',
                    'CB']

#from martini_v2.0_ions.itp
ions_list = ['NA+', 'CL-', 'NC3+', 'CA+']

#from martini_v2.0_sugars.itp
sugars_list =   ['GLUC',
                 'FRUC',
                 'SUCR',
                 'MALT',
                 'CELL',
                 'KOJI',
                 'SOPH',
                 'NIGE',
                 'LAMI',
                 'TREH',
                 'Maltoheptaose']



costituents = {'aminoacids' : set(aminoacids_list),
              'sugars' : set(sugars_list),
              'peptide' :  set(aminoacids_list).union(set(sugars_list))
              }