#! /usr/bin/env python
""" 
    Routines to performe virtual screening
    These routines were developed by:
    Rodrigo Antonio Faccioli - rodrigo.faccioli@usp.br / rodrigo.faccioli@gmail.com  
    Leandro Oliveira Bortot  - leandro.bortot@usp.br / leandro.obt@gmail.com 
"""

import vina
import ConfigParser as configparser
from analysis import call_analysis as ana

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
		print "Analysing "+receptor
		ana.call_vs_analysis()

main()