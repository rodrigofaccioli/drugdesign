import pyspark.sql.functions as func
from pyspark import SparkContext, SparkConf, SparkFiles
from pyspark.sql import SQLContext, Row
import ConfigParser as configparser
import os
from datetime import datetime

def save_result_only_pose(path_file_result_file_only_pose, df_result):
	list_aux = []	
	f_file = open(path_file_result_file_only_pose, "w")
	header = "# poses that were filtered by residues from hydrogen bind\n"
	f_file.write(header)				
	for row in df_result.collect(): 
		if row.pose not in list_aux: #removes duplicates poses
			line = str(row.pose)+"\n"
			f_file.write(line)
			list_aux.append(row.pose)				
	f_file.close()

def save_result(path_file_result_file, df_result):
	f_file = open(path_file_result_file, "w")
	header = "#ligand_atom\taccept_or_donate\treceptor_residue\treceptor_atom\tdistance[A]\tangle[deg]\tpose\n"
	f_file.write(header)				
	for row in df_result.collect():
		line = str(row.ligand_atom)+"\t"+str(row.accept_or_donate)+"\t"+str(row.receptor_residue)+"\t"+str(row.receptor_atom)+"\t"+str(row.distance)+"\t"+str(row.angle)+"\t"+str(row.pose)+"\n"
		f_file.write(line)				
	f_file.close()

def save_log(finish_time, start_time):
	log_file_name = 'hydrogen_bind_residue_selection.log'
	current_path = os.getcwd()
	path_file = os.path.join(current_path, log_file_name)
	log_file = open(path_file, 'w')

	diff_time = finish_time - start_time
	msg = 'Starting ' + str(start_time) +'\n'
	log_file.write(msg)
	msg = 'Finishing ' + str(finish_time) +'\n'
	log_file.write(msg)
	msg = 'Time Execution (seconds): ' + str(diff_time.total_seconds()) +'\n'
	log_file.write(msg)


def main():

	config = configparser.ConfigParser()
	config.read('config.ini')

	#Number of poses to select by buried area
	number_poses_to_select_hydrogen_bind = int(config.get('DRUGDESIGN', 'number_poses_to_select_hydrogen_bind') )
	# list of residues to select buried area
	file_select_hydrogen_bind = config.get('DRUGDESIGN', 'file_residue_to_select_hydrogen_bind')
	#Path that contains all files for analysis
	path_analysis = config.get('DEFAULT', 'path_analysis')	
	#File for saving the filtered buried area
	result_file_to_select_hydrogen_bind = config.get('DRUGDESIGN', 'result_file_to_select_hydrogen_bind')
	#File for saving the filtered buried area only poses
	result_file_to_select_hydrogen_bind_only_pose = config.get('DRUGDESIGN', 'result_file_to_select_hydrogen_bind_only_pose')
	
	# Create SPARK config
	maxResultSize = str(config.get('SPARK', 'maxResultSize'))
	conf = (SparkConf().set("spark.driver.maxResultSize", maxResultSize))

	# Create context
	sc = SparkContext(conf=conf)
	sqlCtx = SQLContext(sc)

	start_time = datetime.now()

	#load all-residue_hbonds_4.0A_30.0deg.dat file
	path_file_hydrogen_bind = os.path.join(path_analysis, "all-residue_hbonds_4.0A_30.0deg.dat")
	all_residue	= sc.textFile(path_file_hydrogen_bind)
	header = all_residue.first() #extract header	

	#Spliting file by \t
	all_residue_split = all_residue.filter(lambda x:x !=header).map(lambda line: line.split("\t"))
	all_residue_split = all_residue_split.map(lambda p: Row( ligand_atom=str(p[0]), accept_or_donate=str(p[1]), receptor_residue=str(p[2]), receptor_atom=str(p[3]), distance=float(p[4]), angle=float(p[5]), pose=str(p[6]) ))

	#Creating all_residue Dataframe
	df_all_residue = sqlCtx.createDataFrame(all_residue_split)	
	df_all_residue.registerTempTable("all_residue")

	#Creating resudue list as Dataframe
	residue_list = sc.textFile(file_select_hydrogen_bind)	
	header = residue_list.first() #extract header		
	#Spliting file by \t
	residue_listRDD = residue_list.filter(lambda x:x !=header).map(lambda line: line)
	residue_listRDD = residue_listRDD.map(lambda p: Row( residue=str(p).strip() ))

	df_residue_list = sqlCtx.createDataFrame(residue_listRDD)	
	df_residue_list.registerTempTable("residue_list")

	#Getting all information based on list of residues
	sql = """
	       SELECT all_residue.*
	       FROM all_residue 
	       JOIN residue_list ON residue_list.residue = all_residue.receptor_residue
	      """
	df_result = sqlCtx.sql(sql)

	#Saving result
	path_file_result_file = os.path.join(path_analysis, result_file_to_select_hydrogen_bind)
	save_result(path_file_result_file, df_result)	

	path_file_result_file_only_pose = os.path.join(path_analysis, result_file_to_select_hydrogen_bind_only_pose)
	save_result_only_pose(path_file_result_file_only_pose, df_result)	

	finish_time = datetime.now()

	save_log(finish_time, start_time)

main()