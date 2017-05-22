# Description
This directory's goal is performing Virtual Screening using AutoDock Vina [1], MGLTools [2], and Apache Spark[3].

All the programs mentioned above were designed for docking, so they had to be expanded to allow for Virtual Screening.
We achieved this by using three different strategies, each of them being stored in a specific directory. 

-   1.mpi
-   2.python
-   3.spark

***************************************************************************************************************************
###1.mpi
* The purpose of mpi is to allow multiple independent instances of Vina to run at the same time. 
See mpi/READMINE for detailed information.

***************************************************************************************************************************
###2.python
* This directory contains Python scripts in order to streamline Vina and provide a brief analysis of the results.
See python/READMINE for detailed information.
 
***************************************************************************************************************************
###3.spark
* Apache Spark [3] is used for more in-depth analysis of the Virtual Screening's results. 
See spark/READMINE for detailed information.
***************************************************************************************************************************


# Installation 


Please download and install the following software:

1) AutoDock Vina - http://vina.scripps.edu/
2) MGLTools-1.5.6 - http://mgltools.scripps.edu/downloads 
3) python-matplotlib package - http://matplotlib.org/users/installing.html
4) cmake - http://www.cmake.org/
5) Massage Passing Interface (MPI) 

Then, compile Virtual Screening:
```
mkdir mpi/build
cd mpi/build
cmake ../
make
```
-----------------------------
##Perform a Virtual Screening

After successful installation, please follow the steps described bellow:

1) Download compounds from ZINC Database: http://zinc.docking.org/
2) Prepare the config.ini file
3) Extract the compounds from ZINC Database
4) Prepare ligands
5) Prepare receptors
6) Prepare AutoDock Vina's Box using Pymol
7) Prepare input files
8) Run the Virtual Screening MPI
9) Evaluate Results

* [1] http://vina.scripps.edu/
* [2] http://mgltools.scripps.edu/
* [3] http://spark.apache.org/

