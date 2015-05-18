#! /usr/bin/env python
""" 
    Routines for ploting energy histogram of virtual screening
    These routines were developed by:
    Rodrigo Antonio Faccioli - rodrigo.faccioli@usp.br / rodrigo.faccioli@gmail.com  
    Leandro Oliveira Bortot  - leandro.bortot@usp.br / leandro.obt@gmail.com 
"""

import sys
import os
import analysis
import histogram_energy as h_eger 
import ConfigParser as configparser


def main():
	config = configparser.ConfigParser()
	config.read('config.ini')

	path_analysis = config.get('DEFAULT', 'path_analysis')	
	h_eger.build_histogram_energy(path_analysis, analysis.get_histogram_filename())	

main()