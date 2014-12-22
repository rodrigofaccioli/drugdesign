#! /usr/bin/env python
""" 
    Routines to performe virtual screening
    These routines were developed by:
    Rodrigo Antonio Faccioli - rodrigo.faccioli@usp.br / rodrigo.faccioli@gmail.com  
    Leandro Oliveira Bortot  - leandro.bortot@usp.br / leandro.obt@gmail.com 
"""

import vina
import ConfigParser as configparser
import vs_preparation as prep
from analysis import call_analysis as ana

def main():
	config = configparser.ConfigParser()
	config.read('config.ini')

	#Preparing the virtual screening enviroment
	prep.preparing_vs(config.get('DEFAULT', 'path_save_structure'),
		config.get('DEFAULT', 'path_save_log'),
		config.get('DEFAULT', 'path_analysis'),
		config.get('DEFAULT', 'pdbqt_ligand_path'),
		config.get('DEFAULT', 'pdbqt_receptor_path')
		)

	#Obtain all receptors to perform the virtual screening
	all_receptor = vina.get_files_pdbqt(config.get('DEFAULT', 'pdbqt_receptor_path'))

	#Performing the the virtual screening for each receptor
	for receptor in all_receptor:
		print receptor
		vina.run_docking(config.get('VINA', 'vina_program'), 
			config.get('VINA', 'config_file'), 
			receptor, 
			config.get('DEFAULT', 'pdbqt_ligand_path'),
			config.get('DEFAULT', 'path_save_structure'),
			config.get('DEFAULT', 'path_save_log')
			)
		print "Analysing "+receptor
		ana.call_vs_analysis()

main()