from pyspark import SparkContext, SparkConf
from pyspark.sql import SQLContext, Row	
import ConfigParser as configparser
import os
from vina_utils import get_file_name_sorted_energy
from datetime import datetime

def save_ligand_efficiency(path_analysis, list_ligand_efficiency):
	path_file_ligand_efficiency = os.path.join(path_analysis, "ligand_efficiency.txt")
	f_ligand_efficiency = open(path_file_ligand_efficiency,"w")
	for ligand_efficiency in list_ligand_efficiency:
		receptor = ligand_efficiency[0]
		ligand = ligand_efficiency[1]
		mode = ligand_efficiency[2]
		lig_eff = "{:10.4f}".format(ligand_efficiency[3])
		line = str(receptor) + "\t" + str(ligand) + "\t" + str(mode) + "\t"+ str(lig_eff) + "\n"
		f_ligand_efficiency.write(line)
	f_ligand_efficiency.close()


def save_ligand_efficiency_log(finish_time, start_time):
	log_file_name = 'vs_ligand_efficiency.log'
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
	
	#Computing ligand efficiency
	ligand_efficiencyRDD = sqlCtx.sql("SELECT vina.receptor, vina.ligand, vina.mode, (vina.energy / database.heavyAtom) as lig_efficiency FROM database JOIN  vina ON vina.ligand = database.ligand ORDER BY lig_efficiency") 
	ligand_efficiencyRDD = ligand_efficiencyRDD.map(lambda p: (p.receptor, p.ligand, p.mode, p.lig_efficiency) ).collect()

	#Saving ligand efficiency file
	save_ligand_efficiency(path_analysis, ligand_efficiencyRDD)

	finish_time = datetime.now()

	save_ligand_efficiency_log(finish_time, start_time)

main()