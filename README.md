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


## Getting started
Is suggested to install the package in a ***conda environment*** using **Anaconda**, due to the *active development status* of *Morphoscanner*.

### Installing Anaconda
The [Anaconda installer](https://www.anaconda.com/distribution/ "Anaconda website") can be downloaded and installed in the user system using [the instructions](https://docs.anaconda.com/anaconda/install/linux/ "Installation Instructions").

### Conda env creation
***Conda envs*** can be created following the [conda documentation](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html "conda envs management").

An *env* called *morphoscanner* can be created with: 

>conda create -n morphoscanner python pip numpy pandas mdanalysis tqdm pytorch networkx

If you have an *Nvidia GPU* you can use hardware acceleration installing the package ***cudatoolkit***.

To date *pytorch* package requires *cudatoolkit=9.2* or *cudatoolkit=10.1*.

You can install the desidered *cudatoolkit* (e.g. *cudatoolkit 10.1*) with:

> conda install cudatoolkit=10.1


You can activate the env with:

> conda activate morphoscanner


### Morphoscanne Installation as Module

Inside the env you have to install morphoscanner:

> pip install git+https://github.com/lillux/sap_analysis.git#egg=morphoscanner

You need to be a collaborator of the project to download the package. The prompt will ask for *username* and *password*.

Then it will be installed in your env. You can now use morphoscanner from your *IDE* or *Python Console*.


### Morphoscanner as Script

You can download morphoscanner from its [Github repository](https://github.com/lillux/sap_analysis "Morphoscanner repository").

It will be downloaded as a compressed archive. Decompress it in a directory. Then move to the directory and launch:

>python main.py

The script will start and the input paths of the trajectory files will be requested.
