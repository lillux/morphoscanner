from setuptools import setup, find_packages

setup(name = 'morphoscanner',
     version = '0.0.2',
     description = 'A library to handle Martini CG Gromacs trajectory data of Self Assembling Peptides',
     url = 'https://github.com/lillux/morphoscanner',
     author = 'Calogero Carlino, Federico Fontana',
     author_email = 'calogero.carlino28@gmail.com',
     license = 'GPLv3',
     zip_safe=False,
     install_requires=['numpy=1.19.1',
                       'scipy=1.5.2',
                       'MDAnalysis=1.0.0',
                       'pandas>=1.0.3',
                       'torch>=1.1.0',
                       'networkx=2.5',
                       'plotly',
                       'matplotlib',
                       'openpyxl',
                       'tqdm'],
                       
      classifiers=[
          'Development Status :: 3 - Alpha',
          'Intended Audience :: Healthcare Industry',
          'Intended Audience :: Education',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
          'Programming Language :: Python :: 3.8',
          'Topic :: Software Development :: Libraries',
          'Topic :: Software Development :: Libraries :: Python Modules',
          'Topic :: Scientific/Engineering :: Bio-Informatics'
      ],
              
      packages=find_packages())

