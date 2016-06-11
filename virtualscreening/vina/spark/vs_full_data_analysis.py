from pyspark import SparkContext, SparkConf
from pyspark.sql import SQLContext, Row	
import ConfigParser as configparser
import os
from vina_utils import get_file_name_sorted_energy
from datetime import datetime

def format_Hydrogen_bind(num_hydrogen_bind):
	if num_hydrogen_bind == None:
		return "0"
	else:
		return str(num_hydrogen_bind)

def save_vs_full_data(path_analysis, list_full_data, full_data_file_name):	
	path_file_vs_full_data = os.path.join(path_analysis, full_data_file_name)
	f_vs_full_data = open(path_file_vs_full_data,"w")
	header = "# affinity[kcal/mol]\tligand_efficiency[kcal/mol.ha]\tnumber_hbonds\tburied_area_lig[nm2]\tburied_area_lig[%]\tburied_area_lig-lig[%]\tburied_area_rec[nm2]\tburied_area_total[nm2]\tpose\n"
	f_vs_full_data.write(header)	
	for full_data in list_full_data:		
		line = str("{:.2f}".format(full_data[0]))+"\t"+str("{:.4f}".format(full_data[1]))+"\t"+format_Hydrogen_bind(full_data[2])+"\t"+str("{:.4f}".format(full_data[3]))+"\t"+str("{:.4f}".format(full_data[4]))+"\t"+str("{:.4f}".format(full_data[5]))+"\t"+str("{:.4f}".format(full_data[6]))+"\t"+str("{:.4f}".format(full_data[7]))+"\t"+ str(full_data[8])+"\n"
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

#**************** Loading file that contains all scores and ligand efficiency
	score_file_name = os.path.join(path_analysis, "summary_energies.dat")
	text_file = sc.textFile(score_file_name)
	header = text_file.first() #extract header		

	#Spliting score file by \t
	rdd_vs_score_sorted_split = text_file.filter(lambda x:x !=header).map(lambda line: line.split("\t"))
	#rdd_vs_score_sorted = rdd_vs_score_sorted_split.map(lambda p: Row(receptor=str(p[0]), ligand=str(p[1]), mode=int(p[2]), energy=float(p[3]) ))
	rdd_vs_score_sorted = rdd_vs_score_sorted_split.map(lambda p: Row(affinity=float(p[0]), ligand_efficiency=float(p[1]), pose=str(p[2]) ))	
	#Creating Vina Datafrase based on score file
	vina_table = sqlCtx.createDataFrame(rdd_vs_score_sorted)	
	vina_table.registerTempTable("vina_lig_efficiency")
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

#**************** Loading Buried Area total
	buried_area_file_name = os.path.join(path_analysis,"summary_buried_area_total.dat")
	buried_area_file = sc.textFile(buried_area_file_name)

	#Spliting file by \t
	header = buried_area_file.first() #extract header		
	rdd_buried_area_split = buried_area_file.filter(lambda x:x !=header).map(lambda line: line.split("\t"))
	#rdd_buried_area = rdd_buried_area_split.map(lambda p: Row( receptor=str(p[0]), ligand=str(p[1]), mode=int(p[2]), buried_lig_rec=float(p[3]), buried_lig_rec_perc=float(p[4]), buried_lig_lig_perc=float(p[5]) ))
	rdd_buried_area = rdd_buried_area_split.map(lambda p: Row( buried_area_total=float(p[0]), pose=str(p[1]) ))

	#Creating buried Dataframe
	buried_table = sqlCtx.createDataFrame(rdd_buried_area)	
	buried_table.registerTempTable("buriedArea_total")
#**************** Finish	

#**************** Loading Buried Area receptor
	buried_area_file_name = os.path.join(path_analysis,"summary_buried_areas_receptor.dat")
	buried_area_file_receptor = sc.textFile(buried_area_file_name)
	header = buried_area_file_receptor.first() #extract header	

	#Spliting file by \t
	buried_area_file_receptor_split = buried_area_file_receptor.filter(lambda x:x !=header).map(lambda line: line.split("\t"))
	buried_area_file_receptor = buried_area_file_receptor_split.map(lambda p: Row( buried_area_receptor=float(p[0]), pose=str(p[1]) ))

	#Creating buried Dataframe
	buried_area_file_receptor_table = sqlCtx.createDataFrame(buried_area_file_receptor)	
	buried_area_file_receptor_table.registerTempTable("buried_area_receptor")
#**************** Finish	

#**************** Loading Buried Area ligand
	buried_area_file_name = os.path.join(path_analysis,"summary_buried_area_ligand.dat")
	buried_area_file_ligand = sc.textFile(buried_area_file_name)
	header = buried_area_file_ligand.first() #extract header	

	#Spliting file by \t
	buried_area_file_ligand_split = buried_area_file_ligand.filter(lambda x:x !=header).map(lambda line: line.split("\t"))
	buried_area_file_ligand = buried_area_file_ligand_split.map(lambda p: Row( buried_area_lig=float(p[0]), buried_area_lig_perc=float(p[1]), buried_area_lig_lig_perc=float(p[2]), pose=str(p[3]) ))

	#Creating buried Dataframe
	buried_area_file_ligand_table = sqlCtx.createDataFrame(buried_area_file_ligand)	
	buried_area_file_ligand_table.registerTempTable("buried_area_ligand")
#**************** Finish	

#**************** Loading Hydrogen Bind 
	hydrogen_bind_num_pose_file_name = os.path.join(path_analysis,"summary_hbonds_4.0_30.0.dat")
	hydrogen_bind_num_pose_file = sc.textFile(hydrogen_bind_num_pose_file_name)
	header = hydrogen_bind_num_pose_file.first() #extract header	

	#Spliting file by \t
	rdd_hydrogen_bind_split = hydrogen_bind_num_pose_file.filter(lambda x:x !=header).map(lambda line: line.split("\t"))
	rdd_hydrogen_bind = rdd_hydrogen_bind_split.map(lambda p: Row(numHydroBind=int(p[0]), pose=str(p[1]) ) )
	#Creating buried Dataframe
	hydrogen_bind_table = sqlCtx.createDataFrame(rdd_hydrogen_bind)	
	hydrogen_bind_table.registerTempTable("hydrogenbind")
#**************** Finish	

	#Creating SQL command
	sql = ""
	sql = "SELECT vina_lig_efficiency.pose, vina_lig_efficiency.affinity, vina_lig_efficiency.ligand_efficiency"
	sql +=" ,buriedArea_total.buried_area_total"
	sql +=" ,buried_area_receptor.buried_area_receptor"
	sql +=" ,buried_area_ligand.buried_area_lig, buried_area_ligand.buried_area_lig_perc, buried_area_ligand.buried_area_lig_lig_perc "
	sql +=" ,hydrogenbind.numHydroBind	"
	sql +=" FROM vina_lig_efficiency"
	sql +=" JOIN buriedArea_total ON buriedArea_total.pose = vina_lig_efficiency.pose"	
	sql +=" JOIN buried_area_receptor ON buried_area_receptor.pose = vina_lig_efficiency.pose"	
	sql +=" JOIN buried_area_ligand ON buried_area_ligand.pose = vina_lig_efficiency.pose"	
	sql +=" LEFT OUTER	"	
	sql +=" JOIN hydrogenbind ON hydrogenbind.pose = vina_lig_efficiency.pose"		
	sql +=" ORDER BY vina_lig_efficiency.pose"

	#Getting all data
	full_dataRDD = sqlCtx.sql(sql) 
	full_dataRDD = full_dataRDD.map(lambda p: (p.affinity, p.ligand_efficiency, p.numHydroBind, p.buried_area_lig, p.buried_area_lig_perc, p.buried_area_lig_lig_perc, p.buried_area_total, p.buried_area_receptor, p.pose) ).collect()

	#Saving file
	save_vs_full_data(path_analysis, full_dataRDD, full_data_file_name)	

	finish_time = datetime.now()

	save_vs_full_data_analysis_log(finish_time, start_time)

main()