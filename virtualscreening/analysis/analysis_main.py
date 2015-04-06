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

	ana.call_vs_analysis(path_analysis, config.get('DEFAULT', 'path_save_log'))

main()