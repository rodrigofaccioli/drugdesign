from pyspark.sql import SQLContext, Row
from vina_utils import get_ligand_from_receptor_ligand_model

"""
	Creates data frame of residue list
	sqlCtx - spark SQL context
	residue_listRDD - RDD for creating data frame. It had been created by load_file_select_hydrogen_bond function	
"""
def create_df_residue_list(sqlCtx, residue_listRDD):
	df_residue_list = sqlCtx.createDataFrame(residue_listRDD)	
	df_residue_list.registerTempTable("residue_list")
	return df_residue_list

"""
	Creates data frame of all residues for hydrogen bond
	sqlCtx - spark SQL context
	residue_listRDD - RDD for creating data frame. It had been created by load_file_all_residue_hbonds function	
"""
def create_df_all_residue(sqlCtx, all_residue_split):
	df_all_residue = sqlCtx.createDataFrame(all_residue_split)	
	df_all_residue.registerTempTable("all_residue")
	return df_all_residue

"""
	Creates data frame of all residues filtered by residue list
	sqlCtx - spark SQL context	
	Important: Before running this function must execute the functions
	create_df_all_residue and create_df_residue_list
"""
def create_df_all_residue_filtered_by_res_list(sqlCtx):
	#Getting all information based on list of residues
	sql = """
	   	SELECT all_residue.*
	   	FROM all_residue 
	   	JOIN residue_list ON residue_list.residue = all_residue.receptor_residue
	    """
	df_result = sqlCtx.sql(sql)
	df_result.registerTempTable("residues_filtered_by_list")
	return df_result

"""
	Group by poses all residues filtered by residue list
	sqlCtx - spark SQL context	
	Important: Before running this function must execute the function
	create_df_all_residue_filtered_by_res_list
"""
def get_group_by_poses_all_residue_filtered_by_res_list(sqlCtx):	
	sql = """
	  	SELECT pose, count(*) as num_res
	   	FROM residues_filtered_by_list 
	    GROUP BY pose
	    ORDER BY num_res DESC 
	   """	
	df_result = sqlCtx.sql(sql)	
	return df_result

"""
	Creates dataframe normalized Hydrogen Bond by donors and acceptors
	sqlCtx - spark SQL context 
    df_only_poses - data frame created by get_group_by_poses_all_residue_filtered_by_res_list function
	Important: 		
	database is created by load_database function from database_io file.
	This load_database function creates RDD only. 
	Therefore, the lines below must be executed before calling this function
    #Loading database
	rdd_database = load_database(sc, ligand_database)
	#Creating Dataframe
	database_table = sqlCtx.createDataFrame(rdd_database)	
	database_table.registerTempTable("database")
"""
def create_df_normalized_by_donors_acceptors(sqlCtx, df_only_poses):	
	normalizedRDD = df_only_poses.map(lambda p: Row(num_res=int(p.num_res), ligand=get_ligand_from_receptor_ligand_model(p.pose), pose=str(p.pose) ) ).collect()
	#Creating Dataframe 
	normalized_residues_filtered_by_list_table = sqlCtx.createDataFrame(normalizedRDD)	
	normalized_residues_filtered_by_list_table.registerTempTable("normalized_residues_filtered_by_list")
	# Normalized Hydrogen Bond by donors and acceptors
	sql = """
		SELECT pose, (b.num_res / a.hb_donors_acceptors) as normalized_hb
		FROM database a 
		JOIN normalized_residues_filtered_by_list b ON b.ligand = a.ligand
		ORDER BY normalized_hb DESC 
	   """
	df_result = sqlCtx.sql(sql)
	return df_result


"""
	Creates dataframe normalized Hydrogen Bond by heavy atoms
	sqlCtx - spark SQL context 
	Important: 		
	database is created by load_database function from database_io file.
	This load_database function creates RDD only. 
	Therefore, the lines below must be executed before calling this function
    #Loading database
	rdd_database = load_database(sc, ligand_database)
	#Creating Dataframe
	database_table = sqlCtx.createDataFrame(rdd_database)	
	database_table.registerTempTable("database")
"""
def create_df_normalized_by_heavy_atoms(sqlCtx):	
	# Normalized Hydrogen Bond by heavy atoms
	sql = """
			SELECT pose, (b.num_res / a.heavyAtom) as normalized_hb
			FROM database a 
			JOIN normalized_residues_filtered_by_list b ON b.ligand = a.ligand
			ORDER BY normalized_hb DESC 
	      """
	df_result = sqlCtx.sql(sql)
	return df_result

"""
	Creates dataframe of hydrogen bond
	sqlCtx - spark SQL context
	rdd_hydrogen_bond - RDD for creating dataframe. It had been created by load_file_summary_hbonds function
"""
def create_df_hydrogen_bond(sqlCtx, rdd_hydrogen_bond):
	hydrogen_bond_table = sqlCtx.createDataFrame(rdd_hydrogen_bond)	
	hydrogen_bond_table.registerTempTable("hydrogenbond")
	return hydrogen_bond_table	
