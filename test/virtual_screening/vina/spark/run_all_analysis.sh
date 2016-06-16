#!/bin/bash

spark_paths=$1
path_files=$2

$spark_paths $path_files"/"create_file_for_analysis.py

$spark_paths $path_files"/"prepare_files_for_analysis.py

$spark_paths $path_files"/"ligand_efficiency.py

$spark_paths $path_files"/"buried_areas.py 0.14 24

$spark_paths $path_files"/"buried_area_receptor.py

$spark_paths $path_files"/"buried_area_ligand.py  0.14 24

$spark_paths $path_files"/"hydrogen_bind.py 4.0 30.0

$spark_paths $path_files"/"vs_full_data_analysis.py