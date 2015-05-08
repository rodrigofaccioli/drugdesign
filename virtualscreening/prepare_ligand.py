#! /usr/bin/env python
""" 
    Routines for preparing ligands: from mol2 to pdbqt files 
    These routines were developed by:
    Rodrigo Antonio Faccioli - rodrigo.faccioli@usp.br / rodrigo.faccioli@gmail.com  
    Leandro Oliveira Bortot  - leandro.bortot@usp.br / leandro.obt@gmail.com 
"""

import ConfigParser as configparser
import vina
import ligand_database as database

def main():
	config = configparser.ConfigParser()
	config.read('config.ini')

	vina.check_for_preparing_ligand(config.get('DEFAULT', 'mol2_path'), 
			config.get('DEFAULT', 'pdbqt_ligand_path'),
			 config.get('VINA', 'pythonsh'),
			 config.get('VINA', 'script_ligand4'))

	vina.prepare_ligand(config.get('DEFAULT', 'mol2_path'), 
			config.get('DEFAULT', 'pdbqt_ligand_path'),
			 config.get('VINA', 'pythonsh'),
			 config.get('VINA', 'script_ligand4'))

	database.build_database(config.get('DEFAULT', 'ligand_database_path_file'), 
		config.get('DEFAULT', 'pdbqt_ligand_path'))	

main()