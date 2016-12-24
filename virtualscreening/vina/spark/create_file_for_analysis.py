from pyspark import SparkContext, SparkConf
import os
import operator
import ConfigParser as configparser
import ntpath
from  vina_utils import get_file_name_sorted_energy, get_files_log, get_separator_filename_mode

def get_log_file_name(myfile):
	""" 
	This function obtains the name of file without filename extension
	"""	
	path, filename = ntpath.split(myfile)
	name =  str(filename.split(".log")[0]) #remove .log
	return name


def build_log_lines(log_file):

	all_log_lines = {}
	log_name = get_log_file_name(log_file)	

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
			#value = (energy, mode)		
			all_log_lines[key] = energy

		if line.find("-----+") >= 0:
			line_for_dict = True
	return all_log_lines

"""
	Create a text file which shows the energies sorted and returns the sorted dictionary
"""	
def create_file_by_sorted_energy(path_analysis, sorted_dic_list):
	text_file = os.path.join(path_analysis,get_file_name_sorted_energy())
	f_file = open(text_file, "w")
	line = "# affinity[kcal/mol]"+"\t"+"pose"+"\n"
	f_file.write(line)
	for l_item in sorted_dic_list:
		#aux = str(l_item[0]).split(get_separator_filename_mode())
		#splited_name = str(aux[0]).split("_-_")
		#receptor = str(splited_name[0])
		#ligand = str(splited_name[1])
		#mode = int(aux[1])
		energy = float(l_item[1])
		#line = str(receptor)+"\t"+str(ligand)+"\t"+str(mode)+"\t"+str(energy)+"\n"
		line = str(energy)+"\t"+str(l_item[0])+"\n"
		f_file.write(line)
	f_file.close()	

"""
	Creates a dictionary from rdd
	all_list_dic is return from collect. It return a list of dictionary
"""	
def create_dictionary_from_rdd(all_list_dic):
	ret_dic = {}
	for list_dic in all_list_dic:
		for key in list_dic:	
			ret_dic[key] = list_dic[key]
	return ret_dic


def main():
	
	sc = SparkContext()

	config = configparser.ConfigParser()
	config.read('config.ini')

	#Broadcast
	path_analysis = config.get('DEFAULT', 'path_analysis')
	path_save_log = config.get('DEFAULT', 'path_save_log')
	path_spark_drugdesign = config.get('DRUGDESIGN', 'path_spark_drugdesign')

	#Adding Python Source file
	sc.addPyFile(os.path.join(path_spark_drugdesign,"vina_utils.py"))
	sc.addPyFile(os.path.join(path_spark_drugdesign,"json_utils.py"))


	#Checking path_analysis
	if not os.path.exists(path_analysis):
		os.makedirs(path_analysis)
	else:
		if len(os.listdir(path_analysis)) > 0:
			raise EnvironmentError("Analysis directory contains files ")

	#preparing log list
	list_obj_log = []
	log_files = get_files_log(path_save_log)
	for flog in log_files:
		list_obj_log.append(flog)

	#appling map and collect
	logRDD = sc.parallelize(list_obj_log)	
	all_lines_dic = logRDD.map(build_log_lines).collect()

	#creating a dictionary from the returned rdd
	dict_from_rdd = create_dictionary_from_rdd(all_lines_dic)
	#sorting dictionary
	sorted_dict_list = sorted(dict_from_rdd.items(), key=operator.itemgetter(1))

	#saving energy file
	create_file_by_sorted_energy(path_analysis, sorted_dict_list)

main()