# Morphoscanner: a Library for the Analysis of Molecular Dynamics Simulations of Self Assembling Peptides

`morphoscanner` is a tool developed to analyze *Gromacs* MD simulations of *SAPs* and recognize specific patterns in the SAPs network, in simulation made with the **Martini CG Force Field**.

`morphoscanner` is able to recognize protein secondary structures, and the emergence of *beta-sheet* structural patterns in systems of peptides. It provides qualitative and quantitative data on the SAP assembling process.

`morphoscanner` is written in Python 3 using an object oriented approach, and it is developed to be versatile and easily accessible.

`morphoscanner` gives to the user a way to handle trajectory data using `torch.tensor`. The software leverages ***parallel computing*** to compute tensor operations. It parallelize operations both on *CPU* and *GPU*, if an **Nvidia GPU** is found on the system and the correct version of *cudatoolkit* is installed.

`morphoscanner` is also an `API`, that enable the user to construct alternative workflow to the one we propose.

**We strongly appreciate feedback from the community, bug reports, and advises on alternative analysis the can be useful to the community. We also strongly appreciate your help in the development!**

## Distribution and intended use

The tool can be distributed using *pip*, *GitHub* repository and used as a ***Python package***.

`morphoscanner` can be imported in an IDE and used to write customized scripts and to perform specific analysis. Morphoscanner can be used to analyze MD trajectory data in a *jupyter-notebook*, and it integrates with the main packages used in the data-science workflow, as Numpy, Pandas, PyTorch, MDAnalysis, Matplotlib and NetworkX.

For a deep review of `morphoscanner` functionalities have a look at the [tutorials](https://github.com/lillux/morphoscanner/blob/main/Tutorials).


## Prerequisites
It is suggested to install the package in a ***conda environment*** using **Anaconda**, due to the *active development status* of *Morphoscanner*.

**If you have an *Nvidia GPU* you can use** *PyTorch* **hardware acceleration by installing the package** *cudatoolkit*.

The *Nvidia Driver*, *cudatoolkit* and *PyTorch* version have to be compatible. The compatibility can be checked in the respective websites:

- [Nvidia Driver and cudatoolkit compatibility](https://docs.nvidia.com/deploy/cuda-compatibility/index.html).
- [PyTorch and cudatoolkit compatibility](https://pytorch.org/).

Tested drivers and packages version are in the following table.

 
 System |Python Version| Nvidia Driver | cudatoolkit | PyTorch
--------|--------------|---------------|-------------|--------
Manjaro 20.1.2 |3.8| 440.100 | 10.2 | 1.6.0
Kubuntu 18.04 |3.7|384.130 | 9.0 | 1.1.0

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

Before creating your `env`, **be sure of the package version** that you need. If you need `pytorch` version `1.1.0`, specify it in the command below, as `pytorch==1.1.0`!!

An *env* called *ms_env* can be created with:
```bash
conda create -n ms_env python=3.8 pip jupyter numpy pandas mdanalysis tqdm pytorch networkx cudatoolkit=10.2 matplotlib scipy plotly
```

The env can be accessed with:
```bash
conda activate ms_env
```

The installed packages can be checked (in the active env) with:
```bash
conda list
```

## Morphoscanner installation

### Morphoscanner Installation as Package

`morphoscanner` has `pygnuplot` package as an optional dependency. If you follow the command below, it **will not** be installed.

Inside the active env, you can install `morphoscanner` with:

```bash
pip install git+https://github.com/lillux/morphoscanner.git#egg=morphoscanner
```

#### Install with optional dependencies
If you want the optional dependency `pygnuplot`, you can install `morphoscanner` with its additional packages using the command below.

```bash
pip install git+https://github.com/lillux/morphoscanner.git#egg=morphoscanner[full]
```

`morphoscanner` will be installed in your env. You can now use `morphoscanner` from your *IDE* or *Python Console*.

#### Install non-default branches
Branches other than the *default branch* can be installed adding the name of the branch that you want to download, like *@branch_name*, after the repository url. For example, to download the ***dev***  branch:

```bash
pip install git+https://github.com/lillux/morphoscanner.git@dev#egg=morphoscanner
```


## Getting started

Using ***Morphoscanner*** as a *Python package* is straightforward, leveraging MDAnalysis capability of I/O.\
The first step is to import `morphoscanner`:
```python
from morphoscanner.trajectory import trajectory
```

The system configuration (.gro in GROMACS) and trajectory files (.xtc or .trr in GROMACS) path is needed:
``` python
_gro = '/path/to/your/gro'
_xtc = '/path/to/your/xtc'
```

Create the *trajectory* class instance:
``` python
trj = trajectory(_gro, _xtc)
```

Multiple consecutive trajectory can be merged and read as a single trajectory:
``` python
trj = trajectory(_gro, (_xtc1, _xtc2, _xtc3))
```

*Specify the  frame sampling*.\
The frame in the trajectory can be sampled.\
To sample all frames just leave `sampling_interval=1`. The value you assign to `sampling_interval` is the number of frame you want to skip for each sampled frame. The value should be an `int`:

``` python
interval = 2
trj.compose_database(sampling_interval = interval)
```

Analyze the simulation dataset (this can take some time):
``` python
trj.analyze_inLoop()
```

Retrieve the data:
``` python
trj.get_data()
```

Show the database with the results of the analysis:
``` python
trj.database
```

A *pandas.DataFrame* will be shown at the end of the analysis.

The database can be saved as an *excel file*, leveraging *pandas* capability:

Set an output path:
```python
output_path = 'path/to/your/directory'
```

Set the name of the file:
```python
file_name = 'name_of_the_output_file'
```

Export the database with .xlsx file extension (you need `openpyxl`):
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

### Complete functionalities
An in deep review of `morphoscanner` functionalities can be found in the [morphoscanner tutorial](https://github.com/lillux/morphoscanner/blob/main/Tutorials/1_morphoscanner_tutorial.ipynb).
