#! /usr/bin/env python

# This script run virtual screening

import vina
import ConfigParser as configparser

def main():
	config = configparser.ConfigParser()
	config.read('config.ini')

	all_receptor = vina.get_files_pdbqt(config.get('DEFAULT', 'pdbqt_receptor_path'))
	for receptor in all_receptor:
		print receptor
	vina.run_docking(config.get('VINA', 'vina_program'), 
		config.get('VINA', 'config_file'), 
		receptor, 
		config.get('DEFAULT', 'pdbqt_ligand_path'),
		config.get('DEFAULT', 'path_save_structure'),
		config.get('DEFAULT', 'path_save_log')
		) 

main()