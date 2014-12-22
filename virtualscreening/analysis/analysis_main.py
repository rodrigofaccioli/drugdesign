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
import histogram_energy as h_eger 
import xvg_histogram_energy_values as xvghist

def main():
	config = configparser.ConfigParser()
	config.read('config.ini') 

	log_dict = analysis.log_files_by_energy(config.get('DEFAULT', 'path_save_log') )

	log_sorted_dict = ana_io.create_file_by_sorted_energy(config.get('DEFAULT', 'path_analysis'), log_dict) 

	xvghist.create_xvg_histogram_energy_values(config.get('DEFAULT', 'path_analysis'), log_sorted_dict)	

	h_eger.build_histogram_energy(config.get('DEFAULT', 'path_analysis'), analysis.get_histogram_filename())

main()