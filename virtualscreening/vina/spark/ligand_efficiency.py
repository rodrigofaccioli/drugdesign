from pyspark import SparkContext, SparkConf
from pyspark.sql import SQLContext, Row	
import ConfigParser as configparser
import os
from vina_utils import get_file_name_sorted_energy, get_ligand_from_receptor_ligand_model, get_name_model_pdb
from datetime import datetime
from database_io import load_database

def save_ligand_efficiency(path_analysis, list_ligand_efficiency):
	path_file_ligand_efficiency = os.path.join(path_analysis, "summary_energies.dat")
	f_ligand_efficiency = open(path_file_ligand_efficiency,"w")
	line = "# affinity[kcal/mol]\tligand_efficiency[kcal/mol.ha]\tpose"+"\n"
	f_ligand_efficiency.write(line)
	for ligand_efficiency in list_ligand_efficiency:
		#receptor = ligand_efficiency[0]
		#ligand = ligand_efficiency[1]
		#mode = ligand_efficiency[2]
		pose = ligand_efficiency[0]
		affinity = ligand_efficiency[1]
		lig_eff = "{:.4f}".format(ligand_efficiency[2])
		#line = str(receptor) + "\t" + str(ligand) + "\t" + str(mode) + "\t"+ str(lig_eff) + "\n"
		line = str(affinity) + "\t" + str(lig_eff) + "\t"+ str(pose) + "\n"
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
	sc.addPyFile(os.path.join(path_spark_drugdesign,"json_utils.py"))

	start_time = datetime.now()

#**************** Loading file that contains all scores
	score_file_name = os.path.join(path_analysis,get_file_name_sorted_energy())
	text_file = sc.textFile(score_file_name)

	#Spliting score file by \t
	header = text_file.first() #extract header
	rdd_vs_score_sorted_split = text_file.filter(lambda x:x !=header)    #filter out header
	rdd_vs_score_sorted_split = rdd_vs_score_sorted_split.map(lambda line: line.split("\t"))
	rdd_vs_score_sorted = rdd_vs_score_sorted_split.map(lambda p: Row(energy=float(p[0]), pose=str(p[1]), ligand=get_ligand_from_receptor_ligand_model(p[1]) )) 
	#Creating Vina Datafrase based on score file
	vina_table = sqlCtx.createDataFrame(rdd_vs_score_sorted)	
	vina_table.registerTempTable("vina")	
#**************** Finish 

#**************** Loading Ligand Database
	rdd_database = load_database(sc, ligand_database)	
	#Creating Dataframe
	database_table = sqlCtx.createDataFrame(rdd_database)	
	database_table.registerTempTable("database")
#**************** Finish 
	
	#Computing ligand efficiency
	ligand_efficiencyRDD = sqlCtx.sql("SELECT vina.pose, vina.energy as affinity, (vina.energy / database.heavyAtom) as lig_efficiency FROM database JOIN  vina ON vina.ligand = database.ligand ORDER BY vina.energy") 
	ligand_efficiencyRDD = ligand_efficiencyRDD.map(lambda p: (p.pose, p.affinity, p.lig_efficiency) ).collect()

	#Saving ligand efficiency file
	save_ligand_efficiency(path_analysis, ligand_efficiencyRDD)

	finish_time = datetime.now()

	save_ligand_efficiency_log(finish_time, start_time)

main()