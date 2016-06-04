#!/bin/bash
gromacs_path=$1
pdb_file=$2
ndx_file=$3
xvg_temp_sasa_lig_pose=$4
probre=$5
ndots=$6
xvg_temp_sasa_lig_complex=$7
pdb_before_vs=$8
xvg_temp_sasa_lig_min=$9

export GMX_MAXBACKUP=-1

# Makes the index file with the ligand (chain z) and the rest (non chain z)
echo -e 'chain z'"\n"'q'"\n" | "$gromacs_path"gmx make_ndx -f $pdb_file -o $ndx_file

# SASA of the isolated ligand in the pose conformation
"$gromacs_path"gmx sasa -f $pdb_file -s $pdb_file -n $ndx_file -o $xvg_temp_sasa_lig_pose -surface "chZ" -output "chZ" -nopbc -noprot -probe $probre -ndots $ndots -xvg "none"

# SASA of the complexed ligand in the pose conformation
"$gromacs_path"gmx sasa -f $pdb_file -s $pdb_file -n $ndx_file -o $xvg_temp_sasa_lig_complex -surface "System" -output "chZ" -nopbc -noprot -probe $probre -ndots $ndots -xvg "none"

# SASA of the isolated ligand in its energy-minimized conformation
"$gromacs_path"gmx sasa -f $pdb_before_vs -s $pdb_before_vs -o $xvg_temp_sasa_lig_min -surface "System" -output "System" -nopbc -noprot -probe $probre -ndots $ndots -xvg "none"

