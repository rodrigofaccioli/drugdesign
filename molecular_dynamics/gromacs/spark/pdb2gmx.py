from pyspark import SparkContext, SparkConf, SparkFiles
from pyspark.sql import SQLContext, Row
from subprocess import Popen, PIPE
from datetime import datetime
import os
import sys
import configparser
from md_description import pdb2gmx
from gromacs_utils import check_file_exists
from os_utils import make_directory, preparing_path, time_execution_log


def load_md_pdb2gmx(file_of_pdb2gmx):
    list_ret = []
    f_file = open(file_of_pdb2gmx, "r")
    count = 0
    for line in f_file:
        splited_line = str(line).split()
        count += 1
        line_number = count
        dir_number = str(line_number)
        pdb_input = str(splited_line[0]).strip()
        obj = pdb2gmx(dir_number, pdb_input)
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

    file_of_pdb2gmx = sys.argv[1]
    check_file_exists(file_of_pdb2gmx)

    start_time = datetime.now()

    # Broadcast
    gromacs_path = sc.broadcast(gromacs_path)
    time_dt = sc.broadcast(time_dt)
    time_dt_pdb = sc.broadcast(time_dt_pdb)

    def run_pdb2gmx(pdb2gmx_obj):

        mddir = pdb2gmx_obj.get_dir_number()
        make_directory(mddir)

        pdb_path = ''.join([os.getcwd(),
                            '/',
                            pdb2gmx_obj.get_pdb_input()])

        work_path = ''.join([os.getcwd(),
                               '/',
                              mddir,
                              '/'])
        os.chdir(work_path)

        command = ''.join(['echo 0 0 0 0 | ',
                           gromacs_path.value,
                           './gmx pdb2gmx ',
                           ' -f ',
                           pdb_path,
                           ' -o ',
                           ' prot ',
                           ' -p ',
                           ' top ',
                           ' -ignh ',
                           ' -ff ',
                           ' amber99sb-ildn ',
                           ' -water ',
                           ' tip3p ',
                           ' -his'])

        proc = Popen(command, shell=True, stdout=PIPE)
        proc.communicate()

    list_obj_pdb2gmx = load_md_pdb2gmx(file_of_pdb2gmx)

    md_pdb2gmxRDD = sc.parallelize(list_obj_pdb2gmx)

    md_pdb2gmxRDD.foreach(run_pdb2gmx)

    finish_time = datetime.now()

    time_execution_log(finish_time, start_time, "gromacs_pdb2gmx.log")

