#! /usr/bin/env python
""" 
    Routines for preparing receptors: from pdb to pdbqt files 
    These routines were developed by:
    Rodrigo Antonio Faccioli - rodrigo.faccioli@usp.br / rodrigo.faccioli@gmail.com  
    Leandro Oliveira Bortot  - leandro.bortot@usp.br / leandro.obt@gmail.com 
"""

import ConfigParser as configparser
import vina


def main():
	config = configparser.ConfigParser()
	config.read('config.ini')

	vina.check_for_preparing_receptor(config.get('DEFAULT', 'pdb_path'), 
			config.get('DEFAULT', 'pdbqt_receptor_path'),
			 config.get('VINA', 'pythonsh'),
			 config.get('VINA', 'script_receptor4'))
	
	vina.prepare_receptor(config.get('DEFAULT', 'pdb_path'), 
			config.get('DEFAULT', 'pdbqt_receptor_path'),
			 config.get('VINA', 'pythonsh'),
			 config.get('VINA', 'script_receptor4'))

main()