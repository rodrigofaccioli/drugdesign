from pyspark.sql import Row	

"""
	load file that contains the residue list for Hydrogen Bond
	sc - spark context
	file_select_hydrogen_bond - path filename of residue list for Hydrogen Bond	
"""
def load_file_select_hydrogen_bond(sc, file_select_hydrogen_bond):
	residue_list = sc.textFile(file_select_hydrogen_bond)	
	header = residue_list.first() #extract header		
	#Spliting file by \t
	residue_listRDD = residue_list.filter(lambda x:x !=header).map(lambda line: line)
	residue_listRDD = residue_listRDD.map(lambda p: Row( residue=str(p).strip() ))
	return residue_listRDD

"""
	load file that contains the residue list for Hydrogen Bond
	sc - spark context
	path_file_hydrogen_bond - path filename of all residues hbond
"""
def load_file_all_residue_hbonds(sc, path_file_hydrogen_bond):	
	all_residue	= sc.textFile(path_file_hydrogen_bond)
	header = all_residue.first() #extract header	

	#Spliting file by \t
	all_residue_split = all_residue.filter(lambda x:x !=header).map(lambda line: line.split("\t"))
	all_residue_split = all_residue_split.map(lambda p: Row( ligand_atom=str(p[0]), accept_or_donate=str(p[1]), receptor_residue=str(p[2]), receptor_atom=str(p[3]), distance=float(p[4]), angle=float(p[5]), pose=str(p[6]) ))
	return all_residue_split

"""
	load file that contains the residue list for Hydrogen Bond
	sc - spark context
	path_file_result_file_only_pose - path filename of only pose file
	Important: In this RDD contains the collumns:
	1) pose - name of pose 
	2) number - number of hydrogen bond 
    3) f_name - file name to save
"""
def load_only_poses_file_hydrogen_bond(sc, path_file_result_file_only_pose):
	only_poseRDD = sc.textFile(path_file_result_file_only_pose)
	header = only_poseRDD.first() #extract header		
	#Spliting file by \t
	only_poseRDD = only_poseRDD.filter(lambda x:x !=header).map(lambda line: line.split("\t"))
	only_poseRDD = only_poseRDD.map(lambda p: Row( pose=str(p[0]).strip(), num_res=int(str(p[1]).strip()), f_name=str(p[1]).strip()+"_hb_"+str(p[0]).strip() ) )
	return only_poseRDD

"""
	load file that contains the residue list for Hydrogen Bond
	sc - spark context
	path_file_result_file_only_pose - path filename of only pose file
	Important: In this RDD contains the collumns:
	1) pose - name of pose 
	2) normalized_hb_res - normalized by number of residue bond 
    3) f_name - file name to save
"""
def load_only_poses_file_hydrogen_bond_normalized_by_residues(sc, path_file_result_file_only_pose):
	only_pose_normalizedRDD = sc.textFile(path_file_result_file_only_pose)
	header = only_pose_normalizedRDD.first() #extract header		
	#Spliting file by \t
	only_pose_normalizedRDD = only_pose_normalizedRDD.filter(lambda x:x !=header).map(lambda line: line.split("\t"))
	only_pose_normalizedRDD = only_pose_normalizedRDD.map(lambda p: Row( pose=str(p[0]).strip(), normalized_hb_res=float(str(p[1]).strip() ), f_name=str(p[1]).strip()+"_hb_"+str(p[0]).strip()  ))
	return only_pose_normalizedRDD	

"""
	load file the summary normalized hydrogen bond files 
	sc - spark context
	path_file_normalized_pose - path filename of normalized pose
"""
def load_file_summary_normalized_hbonds(sc, path_file_normalized_pose):
	normalized_poseRDD = sc.textFile(path_file_normalized_pose)
	header = normalized_poseRDD.first() #extract header		
	#Spliting file by \t
	normalized_poseRDD = normalized_poseRDD.filter(lambda x:x !=header).map(lambda line: line.split("\t"))
	normalized_poseRDD = normalized_poseRDD.map(lambda p: Row( pose=str(p[1]).strip(), normalized=float(str(p[0]).strip()), f_name=str(p[0]).strip()+"_hb_"+str(p[1]).strip() ) )
	return normalized_poseRDD


"""
	load file the summary hydrogen bond file
	sc - spark context
	hydrogen_bond_num_pose_file_name - path filename of hydrogen summary file
"""
def load_file_summary_hbonds(sc, hydrogen_bond_num_pose_file_name):
	hydrogen_bond_num_pose_file = sc.textFile(hydrogen_bond_num_pose_file_name)
	header = hydrogen_bond_num_pose_file.first() #extract header	

	#Spliting file by \t
	rdd_hydrogen_bond_split = hydrogen_bond_num_pose_file.filter(lambda x:x !=header).map(lambda line: line.split("\t"))
	rdd_hydrogen_bond = rdd_hydrogen_bond_split.map(lambda p: Row(numHydroBond=int(p[0]), pose=str(p[1]) ) )
	return rdd_hydrogen_bond

