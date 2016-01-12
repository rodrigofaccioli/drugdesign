#!/bin/bash

#----------------------------------------------------------------------------------
#
#     Script for performance test
#	  Example: ./run_buried_area.sh /home/faccioli/Programs/spark-1.4.1-bin-hadoop2.4/bin/ /home/faccioli/workspace/drugdesign/virtualscreening/spark/ 2 0.14 24
#
#
#                                                         Rodrigo Antonio Faccioli
#                                                                       11 Jan 2016
#-------------------------------------------------------------------------------------

spark_path=$1
drug_design_path=$2
num_core=$3
prob=$4
ndots=$5

i=1
while [ $i -le $num_core ]; do
   "$spark_path"spark-submit --master local[$i] "$drug_design_path""buried_areas.py" $prob $ndots
   mv "vs_buried_areas.log" "vs_buried_areas_"$i".log"   
   let i=$i+1
done	

