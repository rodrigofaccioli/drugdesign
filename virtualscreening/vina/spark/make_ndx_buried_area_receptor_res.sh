gromacs_path=$1
pdb_file=$2
ndx_file_temp_sasa=$3
res=$4

export GMX_MAXBACKUP=-1

echo -e '"rec_'"$res"'" | "lig"'"\n"'q'"\n" | "$gromacs_path"gmx make_ndx -f "$pdb_file" -n $ndx_file_temp_sasa -o $ndx_file_temp_sasa