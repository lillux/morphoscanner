from setuptools import setup, find_packages

setup(name = 'morphoscanner',
     version = '0.0.1',
     description = 'A library to handle gromacs simulation data of Self Assembling Peptides',
     url = 'https://github.com/lillux/sap_analysis',
     author = 'Calogero Carlino, Federico Fontana',
     author_email = 'calogero.carlino28@gmail.com',
     license = 'GPLv3',
     zip_safe=False,
     install_requires=['numpy>=1.18.1',
                       'MDAnalysis>=0.20.1',
                       'pandas>=1.0.1',
                       'pytorch>=1.1.0',
                       'networkx>=2.4',
                       'openpyxl',
                       'tqdm'],
                       
      classifiers=[
          'Development Status :: 3 - Alpha',
          'Intended Audience :: Healthcare Industry',
          'Intended Audience :: Education',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
          'Programming Language :: Python :: 3.7',
          'Topic :: Software Development :: Libraries',
          'Topic :: Software Development :: Libraries :: Python Modules',
          'Topic :: Scientific/Engineering :: Bio-Informatics'
      ],
              
      packages=find_packages())

