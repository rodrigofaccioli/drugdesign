from pyspark import SparkContext, SparkConf, SparkFiles
from pyspark.sql import SQLContext, Row
from subprocess import Popen, PIPE
from datetime import datetime
import os
import sys
from md_description import equilibration
from gromacs_utils import check_file_exists
from os_utils import preparing_path, time_execution_log

try:
    import configparser
except ImportError:
    import ConfigParser as configparser


def load_md_equilibration(file_of_equilibration):
    list_ret = []
    f_file = open(file_of_equilibration, "r")
    for line in f_file:
        splited_line = str(line).split()
        min_all = str(splited_line[0]).strip()
        obj = equilibration(min_all)
        list_ret.append(obj)
    return list_ret

if __name__ == '__main__':
    sc = SparkContext()
    sqlCtx = SQLContext(sc)

    config = configparser.ConfigParser()
    config.read('config.ini')

    # Path for gromacs spark project
    path_spark_md = config.get('DRUGDESIGN', 'path_spark_md')

    # Adding Python Source file
    sc.addPyFile(os.path.join(path_spark_md, "gromacs_utils.py"))
    sc.addPyFile(os.path.join(path_spark_md, "os_utils.py"))
    sc.addPyFile(os.path.join(path_spark_md, "md_description.py"))

    # Path for gromacs program
    gromacs_path = preparing_path(config.get('DRUGDESIGN', 'gromacs_path'))

    time_dt = int(config.get('GROMACS_ANALYSIS', 'time_dt'))
    time_dt_pdb = int(config.get('GROMACS_ANALYSIS', 'time_dt_pdb'))

    file_of_equilibration = sys.argv[1]
    check_file_exists(file_of_equilibration)

    start_time = datetime.now()

    # Broadcast
    gromacs_path = sc.broadcast(gromacs_path)
    time_dt = sc.broadcast(time_dt)
    time_dt_pdb = sc.broadcast(time_dt_pdb)

    def run_equilibration(equilibration_min_obj):

        command = ''.join([gromacs_path.value,
                           './gmx grompp ',
                           ' -f ',
                           ' eq.mdp ',
                           ' -c ',
                           equilibration_min_obj.get_min_all(),
                           ' -p ',
                           'top.top',
                           ' -o ',
                           ' eq.tpr'])
        proc = Popen(command, shell=True, stdout=PIPE)
        proc.communicate()

        command = ''.join(['nohup ',
                           gromacs_path.value,
                          './gmx mdrun ',
                           ' -s ',
                           'eq.tpr',
                           ' -v ',
                           ' -deffnm ',
                           ' eq &'])
        proc = Popen(command, shell=True, stdout=PIPE)
        proc.communicate()

        command = ''.join(['echo 0 | ',
                           gromacs_path.value,
                          './gmx trjconv ',
                           ' -f ',
                           ' eq.gro ',
                           ' -s ',
                           ' eq.tpr ',
                           ' -o ',
                           ' check_eq.pdb ',
                           ' -pbc ',
                           ' mol ',
                           ' -ur ',
                           ' compact'])
        proc = Popen(command, shell=True, stdout=PIPE)
        proc.communicate()

    list_obj_equilibration = load_md_equilibration(file_of_equilibration)

    md_equilibrationRDD = sc.parallelize(list_obj_equilibration)

    md_equilibrationRDD.foreach(run_equilibration)

    finish_time = datetime.now()

    time_execution_log(finish_time, start_time, "equilibration.log")

