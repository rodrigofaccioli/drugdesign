#! /usr/bin/env python
""" 
    Routines for creating the input files that are employed to perform the virtual screening 
    These routines were developed by:
    Rodrigo Antonio Faccioli - rodrigo.faccioli@usp.br / rodrigo.faccioli@gmail.com  
    Leandro Oliveira Bortot  - leandro.bortot@usp.br / leandro.obt@gmail.com 
"""

import os
import vina
import ConfigParser as configparser
import vs_preparation as prep

file_name_docking = "overall_docking_list.txt"
file_name_config  = "config.conf"

def creating_overall_docking_list(current_dir, vina, config):	
	#Obtain all receptors to perform the virtual screening
	all_receptor = vina.get_files_pdbqt(config.get('DEFAULT', 'pdbqt_receptor_path'))
	#Obtain all compounds to perform the virtual screening
	all_compounds = vina.get_files_pdbqt(config.get('DEFAULT', 'pdbqt_ligand_path'))
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
			line = str(only_name_receptor) + "\t" + str(only_name_compound) +"\n"
			file_all_docking.write(line)			
		# File is closed for each receptor because the number of compound can be large. 
		# Therefore, this file is closed to avoid saving amount of data 
		file_all_docking.close()

def valid_end_terminator_path(path):
	if str(path).endswith(os.sep) == False:
		path = str(path)+os.sep
	return path

def creating_config_file(current_dir, config):
	path_file_config = os.path.join(current_dir, file_name_config)
	file_config = open(path_file_config, "w")
	current_dir = valid_end_terminator_path(current_dir)
	line = "Local_Execute = "+ str(current_dir) + "\n"
	file_config.write(line)
	line = "Path_receptor = " + valid_end_terminator_path(str(config.get('DEFAULT', 'pdbqt_receptor_path'))) + "\n"
	file_config.write(line)
	line = "Path_compounds = "+ valid_end_terminator_path(str(config.get('DEFAULT', 'pdbqt_ligand_path'))) + "\n"	
	file_config.write(line)
	line = "Path_out = "+ valid_end_terminator_path(str(config.get('DEFAULT', 'path_save_structure'))) + "\n"	
	file_config.write(line)	
	line = "Path_log = "+ valid_end_terminator_path(str(config.get('DEFAULT', 'path_save_log'))) + "\n"	
	file_config.write(line)		
	path_conf_file = os.path.join(current_dir, config.get('VINA', 'config_file'))
	line = "Config_file = "+ str(path_conf_file) + "\n"	
	file_config.write(line)		
	line = "Vina_program = "+ str(config.get('VINA', 'vina_program')) + "\n"	
	file_config.write(line)	
	line = "Pythonsh = "+ str(config.get('VINA', 'pythonsh')) + "\n"	
	file_config.write(line)	
	line = "Script_ligand4 = "+ str(config.get('VINA', 'script_ligand4')) + "\n"			
	file_config.write(line)	
	line = "Path_mol2 = "+ str(config.get('DEFAULT', 'mol2_path')) + "\n"			
	file_config.write(line)		
	file_config.close()

def main():
	
	config = configparser.ConfigParser()
	config.read('config.ini')

	#Obtaining current directory
	current_dir = os.getcwd()

	#Preparing the virtual screening enviroment
	prep.preparing_vs(config.get('DEFAULT', 'pdbqt_ligand_path'), config.get('DEFAULT', 'pdbqt_receptor_path') )
	vina.check_for_running_docking(config.get('VINA', 'config_file'), config.get('VINA', 'vina_program'))	

	#Creating input files for peforming virtual screening
	creating_overall_docking_list(current_dir, vina, config)
	creating_config_file(current_dir, config)

main()