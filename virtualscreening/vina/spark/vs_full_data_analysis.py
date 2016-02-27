from pyspark import SparkContext, SparkConf
from pyspark.sql import SQLContext, Row	
import ConfigParser as configparser
import os
from vina_utils import get_file_name_sorted_energy
from datetime import datetime

def save_vs_full_data(path_analysis, list_full_data, full_data_file_name):	
	path_file_vs_full_data = os.path.join(path_analysis, full_data_file_name)
	f_vs_full_data = open(path_file_vs_full_data,"w")
	header = ";receptor\tligand\t\tmode\ttorsion\tatom\theavyA\tafinity\tefficiency\tb_lig_rec_perc\tb_lig_lig_perc\n"
	f_vs_full_data.write(header)	
	for full_data in list_full_data:
		line = str(full_data[0])+"\t"+str(full_data[1])+"\t"+str(full_data[2])+"\t"+str(full_data[3])+"\t"+str(full_data[4])+"\t"+str(full_data[5])+"\t"+str(full_data[6])+"\t"+str("{:.4f}".format(full_data[7]))+"\t"+str("{:.4f}".format(full_data[8]))+"\t"+str("{:.4f}".format(full_data[9]))+"\n"
		f_vs_full_data.write(line)
	f_vs_full_data.close()


def save_vs_full_data_analysis_log(finish_time, start_time):
	log_file_name = 'vs_full_data_analysis.log'
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

	sc = SparkContext()
	sqlCtx = SQLContext(sc)

	config = configparser.ConfigParser()
	config.read('config.ini')

	#Path that contains all files for analysis
	path_analysis = config.get('DEFAULT', 'path_analysis')	
	#Ligand Database file
	ligand_database  = config.get('DEFAULT', 'ligand_database_path_file')
	#Path for drugdesign project
	path_spark_drugdesign = config.get('DRUGDESIGN', 'path_spark_drugdesign')

	#Adding Python Source file
	sc.addPyFile(os.path.join(path_spark_drugdesign,"vina_utils.py"))

	#Sufix of completly data file
	full_data_file_name = config.get('DRUGDESIGN', 'full_data_file_name')

	start_time = datetime.now()

#**************** Loading file that contains all scores
	score_file_name = os.path.join(path_analysis,get_file_name_sorted_energy())
	text_file = sc.textFile(score_file_name)

	#Spliting score file by \t
	rdd_vs_score_sorted_split = text_file.map(lambda line: line.split("\t"))
	rdd_vs_score_sorted = rdd_vs_score_sorted_split.map(lambda p: Row(receptor=str(p[0]), ligand=str(p[1]), mode=int(p[2]), energy=float(p[3]) ))
	#Creating Vina Datafrase based on score file
	vina_table = sqlCtx.createDataFrame(rdd_vs_score_sorted)	
	vina_table.registerTempTable("vina")
#**************** Finish 

#**************** Loading Ligand Database
	ligand_database_file = sc.textFile(ligand_database)

	#Spliting score file by \t
	header = ligand_database_file.first() #extract header	
	rdd_database_split = ligand_database_file.filter(lambda x:x !=header).map(lambda line: line.split("\t"))
	rdd_database = rdd_database_split.map(lambda p: Row(ligand=str(p[0]), torsion=int(p[1]), atom=int(p[2]), heavyAtom=int(p[3]) ))
	#Creating Dataframe
	database_table = sqlCtx.createDataFrame(rdd_database)	
	database_table.registerTempTable("database")
#**************** Finish 

#**************** Loading Buried Area file
	buried_area_file_name = os.path.join(path_analysis,"buried_area.txt")
	buried_area_file = sc.textFile(buried_area_file_name)

	#Spliting file by \t
	rdd_buried_area_split = buried_area_file.map(lambda line: line.split("\t"))
	rdd_buried_area = rdd_buried_area_split.map(lambda p: Row( receptor=str(p[0]), ligand=str(p[1]), mode=int(p[2]), sasa_lig_min=float(p[3]), sasa_lig_pose=float(p[4]), sasa_lig_complex=float(p[5]), buried_lig_rec_perc=float(p[6]), buried_lig_lig_perc=float(p[7]) ))
	#Creating buried Dataframe
	buried_table = sqlCtx.createDataFrame(rdd_buried_area)	
	buried_table.registerTempTable("buriedArea")
#**************** Finish	

#**************** Loading Ligand Efficiency
	ligand_efficiency_file_name = os.path.join(path_analysis,"ligand_efficiency.txt")
	ligand_efficiency_file = sc.textFile(ligand_efficiency_file_name)

	#Spliting file by \t
	rdd_ligand_efficiency_split = ligand_efficiency_file.map(lambda line: line.split("\t"))
	rdd_ligand_efficiency = rdd_ligand_efficiency_split.map(lambda p: Row(receptor=str(p[0]),ligand=str(p[1]), mode=int(p[2]), lig_eff=float(p[3]) ) )
	#Creating buried Dataframe
	ligand_efficiency_table = sqlCtx.createDataFrame(rdd_ligand_efficiency)	
	ligand_efficiency_table.registerTempTable("ligandEfficiency")
#**************** Finish	

	#Creating SQL command
	sql = ""
	sql = "SELECT database.ligand, database.torsion, database.atom, database.heavyAtom"
	sql +=" ,vina.receptor, vina.mode, vina.energy"
	sql +=" ,buriedArea.sasa_lig_min, buriedArea.sasa_lig_pose, buriedArea.sasa_lig_complex, buriedArea.buried_lig_rec_perc, buriedArea.buried_lig_lig_perc"
	sql +=" ,ligandEfficiency.lig_eff "
	sql +=" FROM database"
	sql +=" JOIN vina ON vina.ligand = database.ligand"
	sql +=" JOIN buriedArea ON buriedArea.ligand = database.ligand AND "
	sql +="                    buriedArea.receptor = vina.receptor AND"
	sql +="                    buriedArea.mode = vina.mode"
	sql +=" JOIN ligandEfficiency ON ligandEfficiency.ligand = database.ligand AND "
	sql +="                    ligandEfficiency.receptor = vina.receptor AND"
	sql +="                    ligandEfficiency.mode = vina.mode"	
	sql +=" ORDER BY vina.receptor, vina.ligand, vina.mode"

	#Getting all data
	full_dataRDD = sqlCtx.sql(sql) 
	full_dataRDD = full_dataRDD.map(lambda p: (p.receptor, p.ligand, p.mode, p.torsion, p.atom, p.heavyAtom, p.energy, p.lig_eff, p.buried_lig_rec_perc, p.buried_lig_lig_perc) ).collect()

	#Saving file
	save_vs_full_data(path_analysis, full_dataRDD, full_data_file_name)	

	finish_time = datetime.now()

	save_vs_full_data_analysis_log(finish_time, start_time)

main()