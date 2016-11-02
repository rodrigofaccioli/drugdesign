#! /usr/bin/env python
""" 
    Routines to performe virtual screening when MPI is NOT working
    These routines were developed by:
    Rodrigo Antonio Faccioli - rodrigo.faccioli@usp.br / rodrigo.faccioli@gmail.com  
    Leandro Oliveira Bortot  - leandro.bortot@usp.br / leandro.obt@gmail.com 
"""

import os
import vina
import ConfigParser as configparser
import vs_preparation as prep
from analysis import call_analysis as ana

def main():
	config = configparser.ConfigParser()
	config.read((os.getenv('DRUGDESIGN') + '/config/config.ini'))

	#Preparing the virtual screening enviroment
	prep.preparing_vs(config.get('DEFAULT', 'pdbqt_ligand_path'), config.get('DEFAULT', 'pdbqt_receptor_path') )
	vina.check_for_running_docking(config.get('VINA', 'config_file'), config.get('VINA', 'vina_program'))	

	#Obtain all receptors to perform the virtual screening
	all_receptor = vina.get_files_pdbqt(config.get('DEFAULT', 'pdbqt_receptor_path'))

	#Performing the the virtual screening for each receptor
	path_structure_receptor = ""
	path_log_receptor = ""
	path_analysis_receptor = ""
	only_name_receptor = ""
	for receptor in all_receptor:
		print receptor
		#Obtaing paths to the virtual screening enviroment of receptor
		only_name_receptor = os.path.basename(receptor)
		only_name_receptor = str(os.path.splitext(only_name_receptor)[0])	
		path_structure_receptor = os.path.join(config.get('DEFAULT', 'path_save_structure'), only_name_receptor)		
		path_log_receptor = os.path.join(config.get('DEFAULT', 'path_save_log'), only_name_receptor)		
		path_analysis_receptor = os.path.join(config.get('DEFAULT', 'path_analysis'), only_name_receptor)
		#Preparing the virtual screening enviroment of receptor
		prep.preparing_vs_receptor(path_structure_receptor, path_log_receptor, path_analysis_receptor)
		#Running docking receptor against all compounds		
		vina.run_docking(config.get('VINA', 'vina_program'), 
			config.get('VINA', 'config_file'), 
			receptor, 
			config.get('DEFAULT', 'pdbqt_ligand_path'),
			path_structure_receptor,
			path_log_receptor
			)
		print "Analysing "+receptor
		ana.call_vs_analysis(path_analysis_receptor, path_log_receptor)
main()