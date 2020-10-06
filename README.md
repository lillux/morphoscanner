# Morphoscanner: a Library for Self Assembling Peptide Analysis from Molecular Dynamics Simulation made with Gromacs

**Morphoscanner** is a tool developed to analyze *Gromacs* MD simulations of *SAPs* and recognize specific topological patterns in the SAPs network, in simulation made with the **Martini CG Force Field**.

The actual available version of Morphoscanner can recognize *Beta-sheet* topologies and retrieve quantitative and qualitative data on the SAP assembling process.

Morphoscanner is written in Python 3 using an object oriented approach.

Morphoscanner is developed to be versatile and easily accessible at the same time.

The software actually leverages ***parallel computing*** for some steps of the work-flow. It parallelize operations both on *CPU* and *GPU*, if an **Nvidia GPU** is found on the system and the correct version of *cudatoolkit* is installed. Parallelizzation of the full workflow is in development.

The tool can be distributed using *pip* repository and used both as a *script* and as a *python module*.

The script gives ease of usage: it only needs the Gromacs output and some topology info as input, and gives *.xlsx* files (Microsoft Excel) as output.
It can run using the *main.py* file of Morphoscanner.

The python module can be imported in an IDE and used to write customized scripts and to do specific analysis. It can used to analyze MD trajectory data in a jupyter-notebook, and it integrates with the main packages used in the data-science workflow, as numpy, pandas, pytorch, MDAnalysis, matplotlib and NetworkX.


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
pip install git+https://github.com/lillux/morphoscanner.git#egg=morphoscanner
```
You need to be a collaborator of the project to download the package. The prompt will ask for *username* and *password*.

Then it will be installed in your env. You can now use morphoscanner from your *IDE* or *Python Console*.

Other branches other than ***master*** can be installed adding the name of the branch that you want to download, like *@branch_name*, after the repository url:
```bash
pip install git+https://github.com/lillux/morphoscanner.git@branch_name#egg=morphoscanner
```


### Morphoscanner as Script

You can download morphoscanner from its [Github repository](https://github.com/lillux/morphoscanner "Morphoscanner repository").

It will be downloaded as a compressed archive. Decompress it in a directory. Then move to the directory and launch:
```bash
python main.py
```

The script will start and the input paths of the trajectory files will be requested.


## Getting started

Using ***Morphoscanner*** as a *Python module* is straightforward, leveraging MDAnalysis capability of I/O:
``` python
from morphoscanner.trajectory import trajectory
```

The .gro and .xtc or .trr files must be inserted as path:
``` python
_gro = '/path/to/your/gro'
_xtc = '/path/to/your/xtc'
```

Create class instance:
``` python
trj = trajectory(_gro, _xtc)
```

Multiple consecutive trajectory can be merged and read as a single trajectory:
``` python
trj = trajectory(_gro, (_xtc1, _xtc2, _xtc3))
```

*Specify the number of grains composing the peptide backbone.*\

To compose the database the number of aminoacids and other backbone grains has to be known. This information can be obtained in two ways:

- **Automatic Parsing**\
Automatic parsing can be used when the peptides have inserted in the simulation as a separate entities. If this is the case, *Morphoscanner* will parse the .gro file and automatically parse backbone grains and carbohydrate.
To choose *Automatic parsing* just don't insert the *peptide_length* argument, or 

``` python
peptide_length = None
```

- **Assisted Parsing**
If in your simulation you inserted premade aggregate, or groups of peptides as a single entities, you have to specify the size of a single peptide of the group.
To specify the number of grains composing the peptides backbone, just assign the number to the peptide_length argument as an *integer*:

``` python
peptide_length = int
```

*Specify the  frame sampling*.\
The frame in the trajectory can be sampled.\
To sample all frames just leave *interval* argument empty or choose 1 as argument. If you want to sample choose the sampling interval as an *integer*:

``` python
interval = int
```

After choosing the correct *peptide_length* and *interval* (or not, if you want to parse automatically), put them inside the *compose_database* function:
``` python
trj.compose_database(peptide_length = peptide_length, interval = interval)
```

Analyze the database (this can take some time):
``` python
trj.analyze_inLoop()
```

Get data:
``` python
trj.get_data()
```

Get database with results:
``` python
trj.get_database()
```

Show database:
``` python
trj.database
```

A *pandas.DataFrame* will be shown at the end of the analysis.

The database can be saved as an *excel file* leveraging *pandas* capability:

Set an output path:
```python
output_path = 'path/to/your/directory'
```
Set the name of the file:
```python
file_name = 'name_of_the_output_file'
```

Export the database with .xlsx file extension:
```python
trj.database.to_excel(output_path, sheet_name=file_name)
```

### Plotting results

After the analysis the data can be visualized with plotting functions.\
Interactive visualization can be enabled in jupyter-notebook with:
``` python
%matplotlib notebook
```

To deactivate interactive visualization:
``` python
%matplotlib inline
```

- Plot the number of aggregate in the sampled timesteps:
```python
trj.plot_aggregates()
```

- Plot the ratio of contacts antiparallel/(parallel + antiparallel) in the sampled timesteps:
``` python
trj.plot_contacts()
```

- Plot one of the sampled frames, visualizing the aggregate with a color code that define the sense of the majority of the contacts in that aggregate.\
    Green: majority of parallel contacts.\
    Blue: majority of antiparallel contacs.\
    Yellow: equal number of parallel and antiparallel contacts.\
    Gray: no contacts.
```python
trj.plot_frame_aggregate(frame: int)
```
- Plot the graph of one of the sampled frames with qualitative visual indications.\
    Edge thickness: thickness proportional to the number of contacts between the two petides (nodes).\
    Edge green: parallel contact.\
    Edge blue: antiparallel contact.
``` python
trj.plot_graph(0)
```

### Additional data

Additional data can be found in:
``` python
trj.frames[frame]
```

This is a dict that contains a dict for each sampled and analyzed frame, with the data computed during the analysis.