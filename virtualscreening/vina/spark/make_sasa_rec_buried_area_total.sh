#!/bin/bash
gromacs_path=$1
pdb_file=$2
ndx_file=$3
sasa_rec=$4


"$gromacs_path"gmx sasa -f $pdb_file -s $pdb_file -nopbc -n $ndx_file -surface '"System_&_!chZ"' -output '"System_&_!chZ"' -xvg "none" -o $sasa_rec