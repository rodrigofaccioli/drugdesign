from pyspark import SparkContext, SparkConf
from pyspark.sql import SQLContext, Row	
import os

""" This function obtains the name of 
summary statistics  
"""
def get_txt_file_name_summary_statistics():
	return 'vs_energies_summary_statistics.txt'


def get_summary_statistics(sc, rdd_vs_energies_sorted):

	sqlCtx = SQLContext(sc)
	vs_energies_sorted_table = sqlCtx.createDataFrame(rdd_vs_energies_sorted)
	vs_energies_sorted_table.registerTempTable("vs_energies_sorted")

	energies = sqlCtx.sql("SELECT count(energy) as total, min(energy) as min_e, max(energy) as max_e, avg(energy) as avg_e FROM vs_energies_sorted")
	energies_out = energies.map(lambda p: "Total_Model: {},Min_Energy: {},Max_Energy: {},AVG_Energy: {}".format(p.total,p.min_e,p.max_e,p.avg_e)).collect()
	return energies_out

def save_txt_summary_statistics(path_analysis, summary_statistics):
	text_file = os.path.join(path_analysis, get_txt_file_name_summary_statistics())
	f_file = open(text_file, "w")

	for l in summary_statistics[0].split(","):
		line = l+"\n"
		f_file.write(line)
	f_file.close()	

