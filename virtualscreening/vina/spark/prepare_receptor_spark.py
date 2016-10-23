import configparser
import sys
from pyspark import SparkContext, SparkConf, SparkFiles
from pyspark.sql import SQLContext, Row
from datetime import datetime
from os_utils import make_directory, preparing_path, time_execution_log, check_file_exists
from subprocess import Popen, PIPE


def load_receptor_list(file_of_receptor):
    list_ret = []
    f_file = open(file_of_receptor, "r")

    for line in f_file:
        receptor = str(line)
        list_ret.append(receptor)

    return list_ret


if __name__ == '__main__':

    sc = SparkContext()
    sqlCtx = SQLContext(sc)

    config = configparser.ConfigParser()
    config.read('config.ini')

    pythonsh = config.get('VINA', 'pythonsh')
    script_receptor4 = config.get('VINA', 'script_receptor4 ')

    # Broadcast
    pythonsh = sc.broadcast(pythonsh)
    script_receptor4 = sc.broadcast(script_receptor4)

    def run_prepare_receptor_spark(receptor):

        command = ''.join([pythonsh,
                           ' ',
                           script_receptor4,
                           ' -r ',
                           receptor,
                           '.pdb'
                           ' -o ',
                           receptor,
                           '.pdbqt',
                           ' -v '])
        proc = Popen(command, shell=True, stdout=PIPE)
        proc.communicate()


    file_of_receptor = sys.argv[1]
    check_file_exists(file_of_receptor)
    start_time = datetime.now()

    list_receptor = load_receptor_list(file_of_receptor)
    vina_dockingRDD = sc.parallelize(list_receptor)
    vina_dockingRDD.foreach(run_prepare_receptor_spark)

    finish_time = datetime.now()
    time_execution_log(finish_time, start_time, "prepare_receptor_spark.log")

