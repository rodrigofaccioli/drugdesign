#from pyspark import SparkContext, SparkConf
#from pyspark.sql import SQLContext, Row	
import ConfigParser as configparser
import os
from datetime import datetime

""" This function obtains atom number, 3D position 
and atom ID (last column of pdbqt file) from atom list.
"""
def get_lig_values_from_atom_list_2_hydrogen_bind(str_atom_list_from_pdbqt):
	list_return = []
	for atom_line in str_atom_list_from_pdbqt:
		splited_line = atom_line.split()
		atom_num = splited_line[1]
		p_x = splited_line[5]
		p_y = splited_line[6]
		p_z = splited_line[7]
		atomID = str(splited_line[-1])#Get the last column
		item = (atom_num, p_x, p_y, p_z, atomID)
		list_return.append(item)
	return list_return

""" This function obtains atom number, 3D position 
and atom ID (last column of pdbqt file) from atom list.
"""
def get_receptor_values_from_atom_list_2_hydrogen_bind(str_atom_list_from_pdbqt):
	list_return = []
	for atom_line in str_atom_list_from_pdbqt:
		splited_line = atom_line.split()
		res_name = splited_line[3]
		res_num = splited_line[4]
		atom_name = splited_line[2]		
		p_x = splited_line[5]
		p_y = splited_line[6]
		p_z = splited_line[7]
		atomID = str(splited_line[-1])#Get the last column
		item = (res_name, res_num, atom_name, p_x, p_y, p_z, atomID)
		list_return.append(item)
	return list_return


def save_vs_hydrogen_bind_log(finish_time, start_time):
	log_file_name = 'vs_hydrogen_bind.log'
	current_path = os.getcwd()
	path_file = os.path.join(current_path, log_file_name)
	log_file = open(path_file, 'w')

	diff_time = finish_time - start_time
	msg = 'Starting ' + str(start_time) +'\n'
	log_file.write(msg)
	msg = 'Finishing ' + str(finish_time) +'\n'
	log_file.write(msg)
	msg = 'Time Execution (seconds): ' + str(diff_time.total_seconds()) +'\n'
	log_file.write(msg)


def main(): #main()

	sc = SparkContext()
	sqlCtx = SQLContext(sc)

	config = configparser.ConfigParser()
	config.read('config.ini')

	start_time = datetime.now()


	finish_time = datetime.now()

	save_vs_hydrogen_bind_log(finish_time, start_time)

