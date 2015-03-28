drugdesign
==========
came about as an initiative to develop a collaborative project 
in drug design area. 
Nowadays, this project allows to perform virtual screening.
It employs AutoDock Vina and MGLTools-1.5.6 softwares. 
It concedes a virtual screening in which the receptor can be 
both experimental, from PDB database, and theoretical, 
molecular dynamics simulation, for example.


Instalation 
===========
It is necessary to install the follow softwares
1) AutoDock Vina - http://vina.scripps.edu/
2) MGLTools-1.5.6 - http://mgltools.scripps.edu/downloads 
3) python-matplotlib package - http://matplotlib.org/users/installing.html
4) cmake - http://www.cmake.org/
5) Massage Passing Interface (MPI) 
6) Compile virtual screening
6.1) mkdir build
6.2) cd build
6.3) cmake ../
6.4) make

How to run a Virtual Screening?
===============================
After the installing process decribed in Instalation section (see above), 
eight steps are necessary to run virtual screening.

1) Download compounds from ZINC Database: http://zinc.docking.org/
2) Preparing the config.ini file
3) Extract the compounds from ZINC Database
4) Preparing ligands
5) Preparing receptors
6) Preparing box of Auto Vina using Pymol
7) Preparing input files
8) Running the Virtual Screening MPI
9) Analysing the Virtual Screening


Understanding the drugdesign project
====================================
it was decomposed into directory to assist its organization.
In each directory, please, see READMINE file that contains
more information.


