#! /usr/bin/env python

# This script prepares receptors


import ConfigParser as configparser
import vina


def main():
	config = configparser.ConfigParser()
	config.read('config.ini')

	vina.prepare_receptor(config.get('DEFAULT', 'pdb_path'), 
			config.get('DEFAULT', 'pdbqt_receptor_path'),
			 config.get('VINA', 'pythonsh'),
			 config.get('VINA', 'script_receptor4'))

main()