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

def main():

	root_path = sys.argv[1]

	for sub_dir in os.listdir(root_path):
		path_analysis = os.path.join(root_path, sub_dir)
		h_eger.build_histogram_energy(path_analysis, analysis.get_histogram_filename())	

main()