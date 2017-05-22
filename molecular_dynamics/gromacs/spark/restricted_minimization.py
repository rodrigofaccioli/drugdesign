from pyspark import SparkContext, SparkConf, SparkFiles
from pyspark.sql import SQLContext, Row
from subprocess import Popen, PIPE
from datetime import datetime
import os
import sys
import configparser
from md_description import restricted_minimization
from gromacs_utils import check_file_exists
from os_utils import preparing_path, time_execution_log


def load_md_restricted_min(file_of_restricted_min):
    list_ret = []
    f_file = open(file_of_restricted_min, "r")
    count = 0
    for line in f_file:
        splited_line = str(line).split()
        count += 1
        line_number = count
        dir_number = str(line_number)
        min_none = str(splited_line[0]).strip()
        obj = restricted_minimization(dir_number, min_none)
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

    file_of_restricted_min = sys.argv[1]
    check_file_exists(file_of_restricted_min)

    start_time = datetime.now()

    # Broadcast
    gromacs_path = sc.broadcast(gromacs_path)
    time_dt = sc.broadcast(time_dt)
    time_dt_pdb = sc.broadcast(time_dt_pdb)

    def run_restricted_min(restricted_min_obj):

        mddir = restricted_min_obj.get_dir_number()

        min_all_mdp = ''.join([os.getcwd(),
                                '/',
                                'min_all.mdp'])

        work_path = ''.join([os.getcwd(),
                             '/',
                             mddir,
                             '/'])
        os.chdir(work_path)

        command = ''.join([gromacs_path.value,
                           './gmx grompp ',
                           ' -f ',
                           min_all_mdp,
                           ' -c ',
                           restricted_min_obj.get_min_none(),
                           ' -p ',
                           'top.top',
                           ' -o ',
                           ' min_all.tpr'])
        proc = Popen(command, shell=True, stdout=PIPE)
        proc.communicate()

        command = ''.join([gromacs_path.value,
                          './gmx mdrun ',
                           ' -s ',
                           'min_all.tpr',
                           ' -v ',
                           ' -deffnm ',
                           ' min_all'])
        proc = Popen(command, shell=True, stdout=PIPE)
        proc.communicate()

        command = ''.join(['echo 0 | ',
                           gromacs_path.value,
                          './gmx trjconv ',
                           ' -f ',
                           ' min_all.gro ',
                           ' -s ',
                           ' min_all.tpr ',
                           ' -o ',
                           ' check_min_all.pdb ',
                           ' -pbc ',
                           ' mol ',
                           ' -ur ',
                           ' compact'])
        proc = Popen(command, shell=True, stdout=PIPE)
        proc.communicate()

    list_obj_restricted_min = load_md_restricted_min(file_of_restricted_min)

    md_restricted_minRDD = sc.parallelize(list_obj_restricted_min)

    md_restricted_minRDD.foreach(run_restricted_min)

    finish_time = datetime.now()

    time_execution_log(finish_time, start_time, "restricted_minimization.log")

