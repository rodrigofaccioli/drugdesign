# DrugDesign

This project was idealized as a collaborative endeavor in the drug design area where scientists are able to perform drug design pipelines in an easy way. The project was designed to be flexible and scalable, in order to adapt to a variety of different computational environments. The main features until now are described below:

1. Virtual Screening (using one or more receptors).
Analysis can be done by using both experimental (from the PDB database) and theoretical (from molecular dynamics simulation) data, thus providing a broader understanding about the selected compounds and receptors and generating more accurate results.

2. Big Data Analytics through Apache Spark.
Virtual Screening usually outputs a considerable amount of data, by using Spark this data can be analyzed more efficiently.

## Project Structure

The current structure of the project consists in the follow directories:

| Directory | Purpose |
|-----------|---------|
|**`config\`**|Where configuration files are located.|
|**`virtualscreening\`**|Where scripts to perform virtual screening are located. instructions inside of the directory.|
|**`test\`**|Where test scripts are located.|

For more information read instructions inside of the directories.

## Project dependencies

The current version of the project uses python 2.7.x, if the python interpreter isn't built in on your platform you can download it [here](https://www.python.org/downloads/).

## How to install
### Before you start

Some tools are used to automate the installation process. [Fabric](http://www.fabfile.org/) is used to automate the build and deploy process. [Miniconda](http://conda.pydata.org/miniconda.html), a small package manager which is part of [Anaconda](https://www.continuum.io/anaconda-overview), is used to setup the project virtual environment and to install dependencies. The following dependencies are currently used by the project:

1. AutoDock Vina - [vina.scripps.edu](vina.scripps.edu) (Miniconda installation)
2. MGLTools - [mgltools.scripps.edu](mgltools.scripps.edu) (Miniconda installation)
3. Gromacs - [www.gromacs.org](www.gromacs.org) (Miniconda installation)
4. Apache Spark - [spark.apache.org](spark.apache.org) (Fabric installation)
5. OpenMPI - [https://www.open-mpi.org/](https://www.open-mpi.org/) (Fabric installation)

### Installation

To perform the installation process you need only execute the `fabfile.py` script, this will execute the others required steps. If you need more information about how to run this script you can execute the follow command:

`$ python fabfile.py --help`

## How to Contribute

To contribute with project you can see the [issue's](https://github.com/rodrigofaccioli/drugdesign/issues) page and take some stuff to work on, have fun!
