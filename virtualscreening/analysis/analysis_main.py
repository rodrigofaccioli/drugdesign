#! /usr/bin/env python
""" 
    Routines to performe the analysis of virtual screening
    These routines were developed by:
    Rodrigo Antonio Faccioli - rodrigo.faccioli@usp.br / rodrigo.faccioli@gmail.com  
    Leandro Oliveira Bortot  - leandro.bortot@usp.br / leandro.obt@gmail.com 
"""

import ConfigParser as configparser
import analysis 
import analysisio as ana_io

def main():
	config = configparser.ConfigParser()
	config.read('config.ini') 

	log_dict = analysis.log_files_by_energy(config.get('DEFAULT', 'path_save_log') )

	ana_io.create_file_by_sorted_energy(config.get('DEFAULT', 'path_analysis'), 
		log_dict)

main()