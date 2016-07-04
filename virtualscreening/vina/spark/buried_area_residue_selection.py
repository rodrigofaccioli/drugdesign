import pyspark.sql.functions as func
from pyspark import SparkContext, SparkConf, SparkFiles
from pyspark.sql import SQLContext, Row
import ConfigParser as configparser
import os

def save_result_only_pose(path_file_result_file_only_pose, df_result):	
	f_file = open(path_file_result_file_only_pose, "w")
	header = "# poses that were filtered by residues from buried area\n"
	f_file.write(header)				
	for row in df_result.collect():
		line = str(row.pose)+"\n"
		f_file.write(line)				
	f_file.close()

def save_result(path_file_result_file, df_result):
	f_file = open(path_file_result_file, "w")
	header = "# residue\tburied_area_residue[nm2]\tresidue_sasa_buried[%]\tpose\n"
	f_file.write(header)				
	for row in df_result.collect():
		line = str(row.residue)+"\t"+str(row.buried_area_residue)+"\t"+str(row.residue_sasa_buried_perc)+"\t"+str(row.pose)+"\n"
		f_file.write(line)				
	f_file.close()

def main():

	config = configparser.ConfigParser()
	config.read('config.ini')

	#Number of poses to select by buried area
	number_poses_to_select_buried_area = int(config.get('DRUGDESIGN', 'number_poses_to_select_buried_area') )
	# list of residues to select buried area
	file_select_buried_area = config.get('DRUGDESIGN', 'file_residue_to_select_buried_area')
	#Path that contains all files for analysis
	path_analysis = config.get('DEFAULT', 'path_analysis')	
	#File for saving the filtered buried area
	result_file_to_select_buried_area = config.get('DRUGDESIGN', 'result_file_to_select_buried_area')
	#File for saving the filtered buried area only poses
	result_file_to_select_buried_area_only_pose = config.get('DRUGDESIGN', 'result_file_to_select_buried_area_only_pose')
	
	# Create SPARK config
	maxResultSize = str(config.get('SPARK', 'maxResultSize'))
	conf = (SparkConf().set("spark.driver.maxResultSize", maxResultSize))

	# Create context
	sc = SparkContext(conf=conf)
	sqlCtx = SQLContext(sc)

	#load all-residue_buried_areas.dat file
	path_file_buried_area = os.path.join(path_analysis, "all-residue_buried_areas.dat")
	all_residue	= sc.textFile(path_file_buried_area)
	header = all_residue.first() #extract header	

	#Spliting file by \t
	all_residue_split = all_residue.filter(lambda x:x !=header).map(lambda line: line.split("\t"))
	all_residue_split = all_residue_split.map(lambda p: Row( residue=str(p[0]), buried_area_residue=float(p[1]), residue_sasa_buried_perc=float(p[2]), pose=str(p[3]) ))

	#Creating all_residue Dataframe
	df_all_residue = sqlCtx.createDataFrame(all_residue_split)	
	df_all_residue.registerTempTable("all_residue")

	#Creating resudue list as Dataframe
	residue_list = sc.textFile(file_select_buried_area)	
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
	       JOIN residue_list ON residue_list.residue = all_residue.residue
	       ORDER BY all_residue.buried_area_residue DESC
	      """
	df_result = sqlCtx.sql(sql)

	#Saving result
	path_file_result_file = os.path.join(path_analysis, result_file_to_select_buried_area)
	save_result(path_file_result_file, df_result)	


	path_file_result_file_only_pose = os.path.join(path_analysis, result_file_to_select_buried_area_only_pose)
	save_result_only_pose(path_file_result_file_only_pose, df_result)	

main()