""" 
    Routines to create xvg file that represents a histogram of energies from virtual screening execution
    These routines were developed by:
    Rodrigo Antonio Faccioli - rodrigo.faccioli@usp.br / rodrigo.faccioli@gmail.com  
    Leandro Oliveira Bortot  - leandro.bortot@usp.br / leandro.obt@gmail.com 
"""

import os
import analysis

def create_xvg_histogram_energy_values(path_analysis, log_sort_dict):
	"""
    Create a text file which shows the energies sorted and returns the sorted dictionary
    Example:
        >>> create_xvg_histogram_energy_values(path_analysis, log_sort_dict)
    @param path_analysis: place where files are saved
    @type path_analysis: string
    @param log_sort_dict: sorted dictonary of energy values 
    @type log_sort_dict: {}
	"""	
	xvg_file = os.path.join(path_analysis, analysis.get_histogram_filename())	
	dict_file = {}

	ref_energy = float(log_sort_dict[0][1])	
	dict_file[ref_energy] = 0
	for l_item in log_sort_dict:
		if float(l_item[1]) == ref_energy:
			dict_file[ref_energy] = dict_file[ref_energy] + 1
		else:
			ref_energy = float(l_item[1])
			dict_file[ref_energy] = 1

	f_file = open(xvg_file, "w")
	line = "#Energy\tFrequency\n"
	f_file.write(line)
	for key in sorted(dict_file):	
		value = dict_file[key]
		line = str(key)+"\t"+str(value)+"\n"
		f_file.write(line)
	f_file.close()
