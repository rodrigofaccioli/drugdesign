#!/usr/bin/env python
""" 
    Routines for histogram of the virtual screening execution
    These routines were developed by:
    Rodrigo Antonio Faccioli - rodrigo.faccioli@usp.br / rodrigo.faccioli@gmail.com  
    Leandro Oliveira Bortot  - leandro.bortot@usp.br / leandro.obt@gmail.com 
"""

import os
import numpy as np
import matplotlib.pyplot as plt

def build_histogram_energy(path_analysis, xvg_filename):
	"""
    Create a histogram graphic of energies
    Example:
        >>> build_histogram_energy(path_analysis, log_dict)
    @param path_analysis: place where files are saved
    @type path_analysis: string
    @param log_dict: place where log files are
    @type log_dict: {}
	"""	
	x1 = []
	y1 = []	
	xvg_pathfile = os.path.join(path_analysis, xvg_filename)
	f_file = open(xvg_pathfile, "r")

	for line in f_file:
		if str(line).find("#") == -1:
			x1.append( float(str(line).split("\t")[0]) ) 
			y1.append( int(str(line).split("\t")[1]) )			

	fig = plt.figure()

	width = 0.35
	ind = np.arange(len(y1))
	plt.bar(ind, y1)
	plt.xticks(ind + width , x1, rotation=45)

	png_pathfile = str(xvg_pathfile).replace('.xvg','.png')
	plt.savefig(png_pathfile)

