#! /usr/bin/env python
""" 
    Routines to analysis the virtual screening
    These routines were developed by:
    Rodrigo Antonio Faccioli - rodrigo.faccioli@usp.br / rodrigo.faccioli@gmail.com  
    Leandro Oliveira Bortot  - leandro.bortot@usp.br / leandro.obt@gmail.com 
"""

import os
import ConfigParser as configparser
import call_analysis as ana
import docking_time as d_time

def main():
	config = configparser.ConfigParser()
	config.read('config.ini')

	path_analysis = config.get('DEFAULT', 'path_analysis')
	#Checking path_analysis
	if not os.path.exists(path_analysis):
		os.makedirs(path_analysis)
	else:
		if len(os.listdir(path_analysis)) > 0:
			raise EnvironmentError("Analysis directory contains files ")

	#Energy analysis
	ana.call_vs_analysis(path_analysis, config.get('DEFAULT', 'path_save_log'))
	
	#Docking time execution analysis
	path_log_files = os.getcwd()
	d_time.vs_time_docking(path_log_files, path_analysis)

main()