from pyspark import SparkContext, SparkConf, SparkFiles
from pyspark.sql import SQLContext, Row
from subprocess import Popen, PIPE
from datetime import datetime
import os
import sys
import configparser
from md_description import add_ions
from gromacs_utils import check_file_exists
from os_utils import make_directory, preparing_path, time_execution_log


def load_md_add_ions(file_of_add_ions):
    list_ret = []
    f_file = open(file_of_add_ions, "r")
    count = 0
    for line in f_file:
        splited_line = str(line).split()
        count += 1
        line_number = count
        dir_number = str(line_number)
        phys_pattern = str(splited_line[0]).strip()
        obj = add_ions(dir_number, phys_pattern)
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

    file_of_add_ions = sys.argv[1]
    check_file_exists(file_of_add_ions)

    start_time = datetime.now()

    # Broadcast
    gromacs_path = sc.broadcast(gromacs_path)
    time_dt = sc.broadcast(time_dt)
    time_dt_pdb = sc.broadcast(time_dt_pdb)

    def run_add_ions(add_ion_obj):

        mddir = add_ion_obj.get_dir_number()

        min_none_mdp = ''.join([os.getcwd(),
                                '/',
                                'min_none.mdp'])

        work_path = ''.join([os.getcwd(),
                             '/',
                             mddir,
                             '/'])
        os.chdir(work_path)

        command = ''.join([gromacs_path.value,
                           './gmx grompp ',
                           ' -f ',
                           min_none_mdp,
                           ' -c ',
                           ' water_ok.gro ',
                           ' -p ',
                           ' top.top ',
                           ' -o ',
                           'ion_neutral.tpr'])
        proc = Popen(command, shell=True, stdout=PIPE)
        proc.communicate()


        command = ''.join(['echo SOL | ',
                           gromacs_path.value,
                           './gmx genion',
                           ' -s ',
                           ' ion_neutral.tpr ',
                           ' -o ',
                           ' ion_neutral.gro ',
                           ' -p ',
                           ' top.top ',
                           ' -pname ',
                           ' NA ',
                           ' -nname ',
                           ' CL ',
                           ' -pq 1 ',
                           ' -nq -1 ',
                           ' -neutral ',
                           ' -conc ',
                           ' 0.000'])
        proc = Popen(command, shell=True, stdout=PIPE)
        proc.communicate()

        command = ''.join([gromacs_path.value,
                           './gmx grompp ',
                           ' -f ',
                           min_none_mdp,
                           ' -c ',
                           ' ion_neutral.gro ',
                           ' -p ',
                           ' top.top ',
                           ' -o ',
                           'ion_is.tpr'])
        proc = Popen(command, shell=True, stdout=PIPE)
        proc.communicate()

        command = ''.join(['echo SOL | ',
                           gromacs_path.value,
                           './gmx genion ',
                           ' -s ',
                           ' ion_is.tpr ',
                           ' -o ',
                           ' ion_is.gro ',
                           ' -p ',
                           ' top.top ',
                           ' -pname ',
                           ' NA ',
                           ' -nname ',
                           ' CL ',
                           ' -pq 1 ',
                           ' -nq -1 ',
                           ' -neutral ',
                           ' -conc ',
                           add_ion_obj.get_phys_pattern_value()])
        proc = Popen(command, shell=True, stdout=PIPE)
        proc.communicate()

    list_obj_add_ions = load_md_add_ions(file_of_add_ions)

    md_add_ionsRDD = sc.parallelize(list_obj_add_ions)

    md_add_ionsRDD.foreach(run_add_ions)

    finish_time = datetime.now()

    time_execution_log(finish_time, start_time, "gromacs_add_ions.log")

