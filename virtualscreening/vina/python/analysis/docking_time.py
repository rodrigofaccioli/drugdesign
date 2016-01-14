#! /usr/bin/env python
""" 
    Routines to analysis execution time of virtual screening
    These routines were developed by:
    Rodrigo Antonio Faccioli - rodrigo.faccioli@usp.br / rodrigo.faccioli@gmail.com  
    Leandro Oliveira Bortot  - leandro.bortot@usp.br / leandro.obt@gmail.com 
"""

import os
import operator

""" This function obtains all log_docking files 
    in mypath 
"""
def get_files_log_docking(mypath):	
	only_log_file = []
	for root, dirs, files in os.walk(mypath):
		for file in files:
			if file.endswith(".log_docking"):
				f_path = os.path.join(root,file)
				only_log_file.append(f_path)
	return only_log_file


def vs_time_docking(path_log_docking, path_analysis):
	d_docking_time = {}
	str_sep = "__-__"
	file_all_docking_time  = "all_docking_time.txt"
	#Obatain all log_docking files
	all_log_files = get_files_log_docking(path_log_docking)
	#Assign all lines into d_docking_time 
	for log in all_log_files:
		f_file = open(log,"r")
		for line in f_file:			
			line_s = str(line).split("\t")
			receptor = str(line_s[0])
			compound = str(line_s[1])			
			dock_time = float(line_s[2])
			torsion_angle = int(line_s[3])
			num_atom     = int(line_s[4])
			key = receptor+str_sep+compound+str_sep+str(torsion_angle)+str_sep+str(num_atom)
			d_docking_time[key] = dock_time
		f_file.close()
	#Sorting dictionary by dock_time
	sorted_log_dict = sorted(d_docking_time.items(), key=operator.itemgetter(1), reverse=True)

	#Creating file
	f_path = os.path.join(path_analysis,file_all_docking_time)
	f_file = open(f_path,"w")
	f_file.write(";receptor compound torsion_angle atom\n")
	for l_item in sorted_log_dict:
		key_s    = str(l_item[0]).split(str_sep)
		receptor = str(key_s[0])
		compound = str(key_s[1])
		torsion_angle = int(key_s[2])
		num_atom     = int(key_s[3])
		dock_time = float(l_item[1])
		line = str(receptor)+"\t"+str(compound)+"\t"+str(dock_time)+"\t"+str(torsion_angle)+"\t"+str(num_atom)+"\n"
		f_file.write(line)		
	f_file.close()

