import ConfigParser as configparser
import os
import sys
from pyspark import SparkContext, SparkConf, SparkFiles
from pyspark.sql import SQLContext, Row
from datetime import datetime
from os_utils import make_directory, preparing_path, time_execution_log, check_file_exists
from subprocess import Popen, PIPE
from vina_utils import get_files_pdb, get_name_model_pdb


if __name__ == '__main__':

    sc = SparkContext()
    sqlCtx = SQLContext(sc)

    config = configparser.ConfigParser()
    config.read('config.ini')

    pythonsh = config.get('VINA', 'pythonsh')
    script_receptor4 = config.get('VINA', 'script_receptor4')
    pdb_path = config.get('DEFAULT', 'pdb_path')
    pdbqt_receptor_path = config.get('DEFAULT', 'pdbqt_receptor_path')
    path_spark_drugdesign = config.get('DRUGDESIGN', 'path_spark_drugdesign')

    make_directory(pdbqt_receptor_path)

    # Adding Python Source file
    sc.addPyFile(os.path.join(path_spark_drugdesign, "vina_utils.py"))
    sc.addPyFile(os.path.join(path_spark_drugdesign, "json_utils.py"))
    sc.addPyFile(os.path.join(path_spark_drugdesign, "os_utils.py"))

    # Broadcast
    pythonsh = sc.broadcast(pythonsh)
    script_receptor4 = sc.broadcast(script_receptor4)
    pdbqt_receptor_path = sc.broadcast(pdbqt_receptor_path)


    def run_prepare_receptor_spark(receptor):


        receptor_pdbqt = os.path.join(pdbqt_receptor_path.value,
                                      get_name_model_pdb(receptor))

        command = ''.join([pythonsh.value,
                           ' ',
                           script_receptor4.value,
                           ' -r ',
                           receptor,
                           ' -o ',
                           receptor_pdbqt,
                           '.pdbqt',
                           ' -v '])
                                   
        proc = Popen(command, shell=True, stdout=PIPE)
        proc.communicate()


    start_time = datetime.now()

    list_receptor = get_files_pdb(pdb_path)
    vina_dockingRDD = sc.parallelize(list_receptor)
    vina_dockingRDD.foreach(run_prepare_receptor_spark)

    finish_time = datetime.now()
    time_execution_log(finish_time, start_time, "prepare_receptor_spark.log")

