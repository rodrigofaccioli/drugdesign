#!/bin/bash
gromacs_path=$1
pdb_file=$2
ndx_file=$3

echo -e 'chain z'"\n"'q'"\n" | "$gromacs_path"gmx make_ndx -f $pdb_file -o $ndx_file
