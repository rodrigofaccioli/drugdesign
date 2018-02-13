# drugdesign

It was idealized as a collaborative endeavor in the drug design area.

The current build allows users to perform Virtual Screening assays by employing the following software:

1. AutoDock Vina - vina.scripps.edu
2. MGLTools - mgltools.scripps.edu
3. Gromacs - www.gromacs.org
4. Apache Spark - spark.apache.org

The project was designed to be flexible and scalable, in order to adapt to a variety of different computational environments.

The main features are described below:

1. Virtual Screening (using one or more receptors).
Analysis can be done by using both experimental (from the PDB database) and theoretical (from molecular dynamics simulation) data, thus providing a broader understanding about the selected compounds and receptors and generating more accurate results.

2. It selects compounds based on biophysics properties instead of Vina score.

3. Big Data Analytics through Apache Spark.
Virtual Screening usually outputs a considerable amount of data, by using Spark this data can be analyzed more efficiently.
