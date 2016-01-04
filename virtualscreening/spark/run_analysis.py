from pyspark import SparkContext, SparkConf
from pyspark.sql import SQLContext, Row	
import os
import operator
import ConfigParser as configparser
import ntpath
from vina_utils import get_file_name_sorted_energy
from summary_statistics import get_summary_statistics, save_txt_summary_statistics

def main():

	sc = SparkContext()

	config = configparser.ConfigParser()
	config.read('config.ini')

	#Path that contains all files for analysis
	path_analysis = config.get('DEFAULT', 'path_analysis')

	#Path for drugdesign project
	path_spark_drugdesign = config.get('DRUGDESIGN', 'path_spark_drugdesign')

	#Adding Python Source file
	sc.addPyFile(os.path.join(path_spark_drugdesign,"summary_statistics.py"))

	#File that contains sorted energies from all log file
	energy_file_name = os.path.join(path_analysis,get_file_name_sorted_energy())

	text_file = sc.textFile(energy_file_name)

	#Spliting energy file by \t
	rdd_vs_energies_sorted_split = text_file.map(lambda line: line.split("\t"))
	rdd_vs_energies_sorted = rdd_vs_energies_sorted_split.map(lambda p: Row(name=str(p[0]), mode=int(p[1]), energy=float(p[2]) ))

	# Appling Summary and Descriptive Statistics in Energies
	summary_statistics_out = get_summary_statistics(sc, rdd_vs_energies_sorted)
	save_txt_summary_statistics(path_analysis, summary_statistics_out)

main()