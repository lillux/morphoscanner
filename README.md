# Morphoscanner: a Library for Self Assembling Peptide Analysis from Molecular Dynamics Simulation made with Gromacs

**Morphoscanner** is a tool developed to analyze *Gromacs* MD simulations of *SAPs* and recognize specific patterns in the SAPs network, in simulation made with the **Martini CG Force Field**.

The available version of Morphoscanner can recognize *Beta-sheet* topologies and retrieve qualitative and quantitative data on the SAP assembling process.

Morphoscanner is written in Python 3 using an object oriented approach.

Morphoscanner is developed to be versatile and easily accessible at the same time.

The software  leverages ***parallel computing*** to compute tensor operations. It parallelize operations both on *CPU* and *GPU*, if an **Nvidia GPU** is found on the system and the correct version of *cudatoolkit* is installed. Parallelizzation of the full workflow will be added in future developments.

The tool can be distributed using *pip* repository and used as a *python module*.

The script gives ease of usage: it only needs the Gromacs output and some topology info as input, and gives *.xlsx* files (Microsoft Excel) as output.

`morphoscanner` can be imported in an IDE and used to write customized scripts and to do specific analysis. Morphoscanner can used to analyze MD trajectory data in a jupyter-notebook, and it integrates with the main packages used in the data-science workflow, as Numpy, Pandas, PyTorch, MDAnalysis, Matplotlib and NetworkX.


## Prerequisites
Is suggested to install the package in a ***conda environment*** using **Anaconda**, due to the *active development status* of *Morphoscanner*.

**If you have an *Nvidia GPU* you can use** *PyTorch* **hardware acceleration by installing the package** *cudatoolkit*.

The *Nvidia Driver*, *cudatoolkit* and *PyTorch* version have to be compatible. The compatibility can be checked in the respective websites:

- [Nvidia Driver and cudatoolkit compatibility](https://docs.nvidia.com/deploy/cuda-compatibility/index.html).
- [PyTorch and cudatoolkit compatibility](https://pytorch.org/).

Tested drivers and packages version are in the following table.

 
 System | Nvidia Driver | cudatoolkit | PyTorch
--------|---------------|-------------|--------
Manjaro 20.1.2 | 440.100 | 10.2 | 1.6.0
Kubuntu 18.04 | 384.130 | 9.0 | 1.1.0

### Installing Anaconda
The [Anaconda installer](https://www.anaconda.com/distribution/ "Anaconda website") can be downloaded and installed in the user system using [the instructions](https://docs.anaconda.com/anaconda/install/linux/ "Installation Instructions").

### Conda env creation
***Conda envs*** can be created following the [conda documentation](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html "conda envs management").

The channel *conda-forge* is needed to install MDAnalysis.
[Here](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-channels.html) the official *conda* docs on how to manage channels.

Add the conda-forge channel (*--append* will add the channel at the bottom of the channel list, *--add* will add the channel at the top of the channels list).
```bash
conda config --append channels conda-forge
```

An *env* called *morphoscanner* can be created with: 
```bash
conda create -n morphoscanner python=3.8 pip numpy pandas mdanalysis tqdm pytorch networkx cudatoolkit=10.2 matplotlib scipy plotly
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

Branches other than yhe *default branch* can be installed adding the name of the branch that you want to download, like *@branch_name*, after the repository url. For example, to download the v0.0.2 branch:

```bash
pip install git+https://github.com/lillux/morphoscanner.git@v0.0.2#egg=morphoscanner
```


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

*Specify the  frame sampling*.\
The frame in the trajectory can be sampled.\
To sample all frames just leave `sampling_interval=1`. If you want to sample choose the sampling interval as an `int`:

``` python
interval = int
```

Then inside the *compose_database* function:
``` python
trj.compose_database(sampling_interval = interval)
```

Analyze the database (this can take some time):
``` python
trj.analyze_inLoop()
```

Get data:
``` python
trj.get_data()
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

The obtained data can be visualized with plotting functions.\
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
    **Green**: majority of parallel contacts.\
    **Blue**: majority of antiparallel contacs.\
    **Yellow**: equal number of parallel and antiparallel contacts.\
    **Gray**: no contacts.
```python
trj.plot_frame_aggregate(frame: int)
```
- Plot the graph of one of the sampled frames with qualitative visual indications.\
    **Edge thickness**: thickness proportional to the number of contacts between the two petides (nodes).\
    **Edge green**: parallel contact.\
    **Edge blue**: antiparallel contact.
``` python
trj.plot_graph(0)
```

### Additional data

Additional data can be found in:
``` python
trj.frames[frame]
```

This is a dict that contains a dict for each sampled and analyzed frame, with the data computed during the analysis.

### Complete functionality
An in deep review of morphoscanner functionalities can be found in the [`morphoscanner_tutorial.ipynb`]().