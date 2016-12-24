from pyspark import SparkContext, SparkConf, SparkFiles
from pyspark.sql import SQLContext, Row
import os, sys
import ConfigParser as configparser
from database_io import load_database
from database_crud import get_torsion_atom_from_database
from vina_utils import get_files_pdbqt, get_files_pdbqt

file_name_docking = "docking_list.txt"

def creating_docking_list(current_dir, config, sqlCtx):	
	#Obtain all receptors to perform the virtual screening
	all_receptor = get_files_pdbqt( config.get('DEFAULT', 'pdbqt_receptor_path') )
	#Obtain all compounds to perform the virtual screening
	all_compounds = get_files_pdbqt( config.get('DEFAULT', 'pdbqt_ligand_path') )
	#List with all docking
	path_file_docking = os.path.join(current_dir, file_name_docking)	
	file_all_docking = open(path_file_docking, "w")
	line = str(len(all_receptor) * len(all_compounds)) +"\n" #Computes how many docking	will be ran
	file_all_docking.write(line)
	file_all_docking.close()
	for receptor in all_receptor:
		#It is used append mode because, file is closed for each receptor 
		file_all_docking = open(path_file_docking, "a+")
		only_name_receptor = os.path.basename(receptor)
		only_name_receptor = str(os.path.splitext(only_name_receptor)[0])		
		for compound in all_compounds:			
			only_name_compound = os.path.basename(compound)
			only_name_compound = str(os.path.splitext(only_name_compound)[0])			
			torsion_num, atom_num = get_torsion_atom_from_database(sqlCtx, only_name_compound)
			line = str(only_name_receptor) + "\t" + str(only_name_compound) + "\t" + str(torsion_num) + "\t" + str(atom_num) +"\n"
			file_all_docking.write(line)			
		# File is closed for each receptor because the number of compound can be large. 
		# Therefore, this file is closed to avoid saving amount of data 
		file_all_docking.close()

def valid_end_terminator_path(path):
	if str(path).endswith(os.sep) == False:
		path = str(path)+os.sep
	return path


def main():
	
	sc = SparkContext()
	sqlCtx = SQLContext(sc)
	config = configparser.ConfigParser()
	config.read('config.ini')

	#Path where docking list file will be saved
	path_to_save = str(sys.argv[1])

	#Path for drugdesign project
	path_spark_drugdesign = config.get('DRUGDESIGN', 'path_spark_drugdesign')

	sc.addPyFile(os.path.join(path_spark_drugdesign,"database_crud.py"))
	sc.addPyFile(os.path.join(path_spark_drugdesign,"database_io.py"))


#**************** Loading Ligand Database
	ligand_database = config.get('DEFAULT', 'ligand_database_path_file')
	rdd_database = load_database(sc, ligand_database)	
	#Creating Dataframe
	database_table = sqlCtx.createDataFrame(rdd_database)	
	database_table.registerTempTable("database")
#**************** Finish 

	#Creating input files for peforming virtual screening
	creating_docking_list(path_to_save, config, sqlCtx)

main()