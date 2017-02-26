>###Spark directory
In this directory contains Python scripts to execute molecular dynamics and others analysis by Apache Spark [1].

>#####Running molecular dynamic
* pdb2gmx.py (required). Generates .gro and .top

* genbox.py (required). Generate the simulation box for PBC use.

* add_ions.py (required). 
Adds ions by doing the steps of system neutralization and reaching ionic strength required.

* unrestricted_minimization.py(required). Performs unrestricted minimization. 
This step requires clustering.

* restricted_minimization.py (required). Performs restricted minimization. Warning: The output from the previous step is a prerequisite for this step. This step requires clustering.

* equilibration.py (required).

* molecular_dynamic.py (required).

>#####Running analysis

* trajectory.py

* protein_flexibility.py