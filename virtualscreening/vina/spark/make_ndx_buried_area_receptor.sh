#!/bin/bash
gromacs_path=$1
pdb_file=$2
ndx_file_temp_index_z=$3
ndx_file_temp=$4

export GMX_MAXBACKUP=-1

# Creates a selection with chain Z (which contains the ligand). This was done because I couldn't make "chain Z" work with "$path_gmx""/"./gmx select
echo -e 'chain Z'"\n"'"System" & ! "chZ"'"\n"'q'"\n" | "$gromacs_path"gmx make_ndx -f $pdb_file -o $ndx_file_temp_index_z

# Creates a selection with the residues that are closer than 6A to the ligand
"$gromacs_path"gmx select -f "$pdb_file"  -s "$pdb_file" -n "$ndx_file_temp_index_z" -select ' "contact" same residue as (group "System_&_!chZ" and within 0.6 of group "chZ") or (group "chZ")' -on $ndx_file_temp

echo -e 'q'"\n" | "$gromacs_path"gmx make_ndx -f "$pdb_file" -n $ndx_file_temp -o $ndx_file_temp

echo -e 'name 0 complex'"\n"'0 & ! chain Z'"\n"'name 1 rec'"\n"'chain Z'"\n"'name 2 lig'"\n"'splitres 1'"\n"'q' | "$gromacs_path"gmx make_ndx -f "$pdb_file" -o $ndx_file_temp -n $ndx_file_temp
