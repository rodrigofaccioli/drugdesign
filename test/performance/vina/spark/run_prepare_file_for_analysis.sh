#!/bin/bash

#----------------------------------------------------------------------------------
#
#     Script for performance test
#	  Example: ./run_prepare_file_for_analysis.sh /home/faccioli/Programs/spark-1.4.1-bin-hadoop2.4/bin/ /home/faccioli/workspace/drugdesign/virtualscreening/spark/ 2
#
#
#                                                         Rodrigo Antonio Faccioli
#                                                                       11 Jan 2016
#-------------------------------------------------------------------------------------

spark_path=$1
drug_design_path=$2
num_core=$3

i=1
while [ $i -le $num_core ]; do
   "$spark_path"spark-submit --master local[$i] "$drug_design_path""prepare_files_for_analysis.py" 
   mv "vs_prepare_files_for_analysis.log" "vs_prepare_files_for_analysis_"$i".log"
   rm -rf analysis/pdb*
   let i=$i+1
done	
