#!/bin/bash

spark_paths=$1
path_files=$2
d_memory=$3"g"

$spark_paths --driver-memory $d_memory $path_files"/"create_file_for_analysis.py

$spark_paths --driver-memory $d_memory $path_files"/"prepare_files_for_analysis.py

$spark_paths --driver-memory $d_memory $path_files"/"ligand_efficiency.py

$spark_paths --driver-memory $d_memory $path_files"/"buried_areas.py 0.14 24

$spark_paths --driver-memory $d_memory $path_files"/"buried_area_receptor.py

$spark_paths --driver-memory $d_memory $path_files"/"buried_area_ligand.py  0.14 24

$spark_paths --driver-memory $d_memory $path_files"/"hydrogen_bond.py 4.0 30.0

$spark_paths --driver-memory $d_memory $path_files"/"vs_full_data_analysis.py

$spark_paths --driver-memory $d_memory $path_files"/"buried_area_residue_selection.py

$spark_paths --driver-memory $d_memory $path_files"/"hydrogen_bond_residue_selection.py

$spark_paths --driver-memory $d_memory $path_files"/"mult_objective_selection.py
