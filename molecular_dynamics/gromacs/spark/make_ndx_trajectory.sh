#!/bin/bash
gromacs_path=$1
tpr_file=$2
ndx_file=$3

echo -e 'q'"\n" | "$gromacs_path"gmx make_ndx -f $tpr_file -o $ndx_file >/dev/null 2>/dev/null