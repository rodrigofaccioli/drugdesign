""" 
    Routines for analysing the virtual screening execution
    These routines were developed by:
    Rodrigo Antonio Faccioli - rodrigo.faccioli@usp.br / rodrigo.faccioli@gmail.com  
    Leandro Oliveira Bortot  - leandro.bortot@usp.br / leandro.obt@gmail.com 
"""

import os
import ntpath
import operator


def get_separator_filename_mode():
	return '+----+'

""" This function obtains all log files 
    in mypath 
"""
def get_files_log(mypath):
	only_log_file = []
	for root, dirs, files in os.walk(mypath):
		for f in files:
			if f.endswith(".log"):
				f_path = os.path.join(root,f)
				only_log_file.append(f_path)			
	return only_log_file

""" This function obtains the name of file  
	without filename extension.
"""
def get_log_file_name(myfile):
	path, filename = ntpath.split(myfile)
	name =  str(filename.split(".")[0]) #remove .log
	return name

def set_energies_dic_return(log_file, log_name, dic_return):
	"""
    Set dictionary from all log files
    Example:
        >>> set_energies_dic_return(log_file, log_name, dic_return)
    @param log_file: full name of log file
    @type log_file: string
    @param log_name: name of log file without filename extension and path. 
    It is used for dictionary's key
	@type log_name: string
    @param dic_return: dictionary that is stored energies of log file
	@type dic_return: {}    
	"""	
	line_for_dict = False
	f_log = open(log_file, "r")
	for line in f_log:
		if line.find("Writing output") >= 0:
			line_for_dict = False

		if 	line_for_dict == True:
			splited_line = line.split()				
			mode = splited_line[0]
			energy = float(splited_line[1])
			key = log_name+get_separator_filename_mode()+mode
			dic_return[key] = energy

		if line.find("-----+") >= 0:
			line_for_dict = True

def log_files_by_energy(mylog):
	"""
    Create a dictionary from all log files.
    Example:
        >>> ret_dict=log_files_by_energy(mylog)        
    @param mylog: place where log files are
    @type mylog: string
    @return: a dictionary
    @rtype: {}
	"""	
	dic_return = {}
	files = get_files_log(mylog)
	for f in files:
		log_name = get_log_file_name(f)
		set_energies_dic_return(f, log_name, dic_return)		
	return dic_return

