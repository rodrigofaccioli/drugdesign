""" 
    Routines to create sorted energy file from virtual screening execution
    These routines were developed by:
    Rodrigo Antonio Faccioli - rodrigo.faccioli@usp.br / rodrigo.faccioli@gmail.com  
    Leandro Oliveira Bortot  - leandro.bortot@usp.br / leandro.obt@gmail.com 
"""

import os
import operator
import analysis as ana


def create_file_by_sorted_energy(path_analysis, log_dict):
	"""
    Create a text file which shows the energies sorted and returns the sorted dictionary
    Example:
        >>> sorted_dict = create_file_by_sorted_energy(path_analysis, log_dict)
    @param path_analysis: place where files are saved
    @type path_analysis: string
    @param log_dict: dictionary with energy values
    @type log_dict: {}
    @return: sorted dictonary of energy values 
    @rtype: {}
	"""	
	text_file = os.path.join(path_analysis,'vs_energies_sorted.txt')
	f_file = open(text_file, "w")
	f_file.write('Name'+"\t"+'Mode'+"\t"+'Energy'+"\n")

	sorted_log_dict = sorted(log_dict.items(), key=operator.itemgetter(1))
	 
	n = 0
	for l_item in sorted_log_dict:
		aux = str(sorted_log_dict[n][0]).split(ana.get_separator_filename_mode())
		name = str(aux[0])		
		mode = int(aux[1])
		energy = float(sorted_log_dict[n][1])		
		line = str(name)+"\t"+str(mode)+"\t"+str(energy)+"\n"
		f_file.write(line)
		n = n + 1
	f_file.close()	

	return sorted_log_dict