#! /usr/bin/env python

# This script prepares ligands


import ConfigParser as configparser
import vina


def main():
	config = configparser.ConfigParser()
	config.read('config.ini')

	vina.prepare_ligand(config.get('DEFAULT', 'mol2_path'), 
			config.get('DEFAULT', 'pdbqt_ligand_path'),
			 config.get('VINA', 'pythonsh'),
			 config.get('VINA', 'script_ligand4'))


main()