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

[Fabric](http://www.fabfile.org/) is used to automate the build and deploy process. The following dependencies are currently used by the project:

1. AutoDock Vina - [vina.scripps.edu](vina.scripps.edu)
2. MGLTools - [mgltools.scripps.edu](mgltools.scripps.edu)
3. Gromacs - [www.gromacs.org](www.gromacs.org)
4. Apache Spark - [spark.apache.org](spark.apache.org)
5. OpenMPI - [https://www.open-mpi.org/](https://www.open-mpi.org/)

### Installation

To perform the installation process, local or remote, you need only execute the following command:

`$ fab -H hosts -u username -p password --sudo-password=sudo_password build_drugdesign`

For more information about the fab program read the [docs](http://docs.fabfile.org/en/1.12/).

### Execution

To execute drugdesign project you'll need to run `load_variables.sh` before:

`$ ./load_variables.sh`

## How to Contribute

To contribute with project you can see the [issue's](https://github.com/rodrigofaccioli/drugdesign/issues) page and take some stuff to work on, have fun!
