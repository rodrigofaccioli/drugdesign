#!/bin/bash
gromacs_path=$1
xtc_file=$2
tpr_file=$3
ndx_file=$4
select_string=$5
time_dt=$6

"$gromacs_path"gmx select -f "$xtc_file" -s "$tpr_file" -on "$ndx_file" -select "$select_string" -dt $time_dt  >/dev/null 2>/dev/null
