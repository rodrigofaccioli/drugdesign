from pyspark import SparkContext, SparkConf, SparkFiles
from pyspark.sql import SQLContext, Row
import configparser
import os
from vina_utils import get_files_log
from json_utils import create_jsondata_from_docking_output_file, create_json_file
from os_utils import make_directory


def log_to_json(log_file):     
    
    path_to_save = path_to_save_b.value

    log_file_name = str(log_file)
    splited = log_file_name.split('/')
    json_name = splited[-1].replace('log', 'json')
    key_name = json_name.split('.')[0]
    json_data = create_jsondata_from_docking_output_file(log_file)
    json_final_data = {key_name: json_data}
    json_name = os.path.join(path_to_save, json_name )
    create_json_file(json_name, json_final_data)

if __name__ == '__main__':

    config = configparser.ConfigParser()
    config.read('config.ini')

    sc = SparkContext()
    sqlCtx = SQLContext(sc)

    log_file = config.get('DEFAULT', 'path_save_log')    
    path_to_save = config.get('DEFAULT', 'json_log')    

    path_spark_drugdesign = config.get('DRUGDESIGN', 'path_spark_drugdesign')   
    sc.addPyFile(os.path.join(path_spark_drugdesign,"vina_utils.py"))
    sc.addPyFile(os.path.join(path_spark_drugdesign,"json_utils.py"))
    sc.addPyFile(os.path.join(path_spark_drugdesign,"os_utils.py"))    

    #Broadcast
    log_file_b = sc.broadcast(log_file)
    path_to_save_b = sc.broadcast(path_to_save)    

    make_directory(path_to_save)
    all_log_files = get_files_log(log_file)
    log_filesRDD = sc.parallelize(all_log_files)
    log_filesRDD.foreach(log_to_json)
