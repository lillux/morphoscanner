# Morphoscanner: a Library for Self Assembling Peptide Analysis from Molecular Dynamics Simulation made with Gromacs

**Morphoscanner** is a tool developed to analyze *Gromacs* MD simulations of *SAPs* and recognize specific topological patterns in the SAPs network.

The actual available version of Morphoscanner can recognize *Beta-sheet* topologies and retrieve quantitative and qualitative data on the SAP assembling process.

Morphoscanner is written in Python 3 using an object oriented approach.

Morphoscanner is developed to be versatile and easily accessible at the same time.

The software actually leverages ***parallel computing*** for some steps of the work-flow. It parallelize operations both on *CPU* and *GPU*, if an **Nvidia GPU** is found on the system and the correct version of *cudatoolkit* is installed. Parallelizzation of the full workflow is in development.

The tool can be distributed using *pip* repository and used both as a *script* and as a *python module*.

The script gives ease of usage: it only needs the Gromacs output and some topology info as input, and gives *.xlsx* files (Microsoft Excel) as output.
It can be runned using the *main.py* file of Morphoscanner.

The python module can be imported in an IDE and used to write customized scripts and to do specific analysis. It can used to analyze MD trajectory data in a jupyter-notebook, and it integrates with the main packages used in a data-science, as numpy, pandas, pytorch, MDAnalysis, matplotlib and NetworkX.


## Prerequisites
Is suggested to install the package in a ***conda environment*** using **Anaconda**, due to the *active development status* of *Morphoscanner*.

### Installing Anaconda
The [Anaconda installer](https://www.anaconda.com/distribution/ "Anaconda website") can be downloaded and installed in the user system using [the instructions](https://docs.anaconda.com/anaconda/install/linux/ "Installation Instructions").

### Conda env creation
***Conda envs*** can be created following the [conda documentation](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html "conda envs management").

**If you have an *Nvidia GPU* you can use** *pytorch* **hardware acceleration installing the package** *cudatoolkit*.\
**To date** *pytorch* **package requires** *cudatoolkit=9.2* or *cudatoolkit=10.1*.\
**Specify the** *cudatoolkit* **that you want to use in the command below**.\
*cudatoolkit 10.1* **is tested and working**.

An *env* called *morphoscanner* can be created with: 
```bash
conda create -n morphoscanner python pip numpy pandas mdanalysis tqdm pytorch networkx cudatoolkit=10.1
```

The installed packages can be checked (in the active env) with:
```bash
conda list
```

You can activate the env with:
```bash
conda activate morphoscanner
```

## Morphoscanner installation

### Morphoscanner Installation as Module

Inside the env you have to install morphoscanner:

```bash
pip install git+https://github.com/lillux/sap_analysis.git#egg=morphoscanner
```
You need to be a collaborator of the project to download the package. The prompt will ask for *username* and *password*.

Then it will be installed in your env. You can now use morphoscanner from your *IDE* or *Python Console*.


### Morphoscanner as Script

You can download morphoscanner from its [Github repository](https://github.com/lillux/sap_analysis "Morphoscanner repository").

It will be downloaded as a compressed archive. Decompress it in a directory. Then move to the directory and launch:

```bash
python main.py
```

The script will start and the input paths of the trajectory files will be requested.


## Getting started

Using ***Morphoscanner*** as a *Python module* is straightforward:

``` python
from morphoscanner.trajectory import trajectory

#The .gro and .xtc files must be inserted as path:
_gro = '/path/to/your/gro'
_xtc = '/path/to/your/xtc'

#Create class instance
trj = trajectory(_gro, _xtc)


#Compose the database
    # peptide_length (integer) is the number of aminoacid in a peptide.
        # leave it blank to follow topology
        # insert an integer if your want to specify the lenght of the peptide 
        # (all the peptides are of the same length)
        # the integer is the number of aminoacids of one peptide
        # if all the peptides are of the same length.
    #interval (integer) is the interval between sampled frame
        # (default) leave it blank to sample each frame
        # Insert an integer

peptide_length = integer
interval = integer

trj.compose_database(peptide_length = peptide_length, interval = interval)

#Analyze the database.
trj.analyze_inLoop()

#Get data
trj.get_data()

#Get database with results:
trj.get_database()

#Show database
trj.database

```

A *pandas.DataFrame* will be shown at the end of the analysis.

The database can be saved as an *excel file* leveraging *pandas* capability:

```python

# set an output path
output_path = 'path/to/your/directory'

# set the name of the file
file_name = 'name_of_the_output_file'

# export the database with .xlsx file extension
trj.database.to_excel(output_path, sheet_name=file_name)


```