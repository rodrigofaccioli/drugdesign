from pyspark import SparkContext, SparkConf
import ConfigParser as configparser
from subprocess import Popen, PIPE
import os
import vina_utils

def get_name_ligand_pdbqt(reference):	
	name =  os.path.basename(reference)
	if str(name).find(".pdb") > 0:
		name = str(name).replace(".pdb", ".pdbqt")	
	if str(name).find(".mol2") > 0:
		name = str(name).replace(".mol2", ".pdbqt")	
	return name	


""" This function converts mol2 or pdb files 
    to pdbqt files for ligand
"""
def prepare_ligand(lig_vina):
		fpdbqt_filename = get_name_ligand_pdbqt(lig_vina[3])
		fpdbqt = os.path.join(lig_vina[0],fpdbqt_filename)				
		process = Popen([lig_vina[1], lig_vina[2], '-l', lig_vina[3], '-v', '-o', fpdbqt, '-U', 'nphs_lps', '-A', 'hydrogens'], stdout=PIPE, stderr=PIPE)		
		stdout, stderr = process.communicate()	 	

def save_log(finish_time, start_time):
	log_file_name = 'prepare_ligand.log'
	current_path = os.getcwd()
	path_file = os.path.join(current_path, log_file_name)
	log_file = open(path_file, 'w')

	diff_time = finish_time - start_time
	msg = 'Starting ' + str(start_time) +'\n'
	log_file.write(msg)
	msg = 'Finishing ' + str(finish_time) +'\n'
	log_file.write(msg)
	msg = 'Total ' + str(diff_time) +'\n'
	log_file.write(msg)


def main():
	
	sc = SparkContext()
	config = configparser.ConfigParser()
	config.read('config.ini')

	#Broadcast - global
	path_pdbqt     = config.get('DEFAULT', 'pdbqt_ligand_path')
	pythonsh       = config.get('VINA', 'pythonsh')
	script_ligand4 = config.get('VINA', 'script_ligand4')

	if not os.path.isdir(path_pdbqt):
		os.mkdir(path_pdbqt)

	start_time = vina_utils.get_time() 

	#preparando lista
	list_obj_lig_vina = []
	mol2_files = vina_utils.get_files_mol2(config.get('DEFAULT', 'mol2_path'))
	for fmol2 in mol2_files:
		obj_lig_vina = (path_pdbqt, pythonsh,script_ligand4, fmol2)
		list_obj_lig_vina.append(obj_lig_vina)

	molRDD = sc.parallelize(list_obj_lig_vina)	
	molRDD.foreach(prepare_ligand)

	finish_time = vina_utils.get_time()

	save_log(finish_time, start_time)

main()




