from pyspark import SparkContext, SparkConf, SparkFiles
from pyspark.sql import SQLContext, Row
from subprocess import Popen, PIPE
from datetime import datetime
import os
import sys
import configparser
from md_description import genbox
from gromacs_utils import check_file_exists
from os_utils import make_directory, preparing_path, time_execution_log


def load_md_genbox(file_of_genbox):
    list_ret = []
    f_file = open(file_of_genbox, "r")
    count = 0
    for line in f_file:
        splited_line = str(line).split()
        count += 1
        line_number = count
        dir_number = str(line_number)
        box_size = str(splited_line[0]).strip()
        obj = genbox(dir_number, box_size)
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

    file_of_genbox = sys.argv[1]
    check_file_exists(file_of_genbox)

    start_time = datetime.now()

    # Broadcast
    gromacs_path = sc.broadcast(gromacs_path)
    time_dt = sc.broadcast(time_dt)
    time_dt_pdb = sc.broadcast(time_dt_pdb)

    def run_genbox(genbox_obj):

        mddir = genbox_obj.get_dir_number()

        min_none_mdp = ''.join([os.getcwd(),
                                '/',
                                'min_none.mdp'])

        work_path = ''.join([os.getcwd(),
                             '/',
                             mddir,
                             '/'])
        os.chdir(work_path)

        command = ''.join(['echo Protein | ',
                           gromacs_path.value,
                          './gmx editconf',
                           ' -f ',
                           'prot.gro',
                           ' -o ',
                           ' box.gro ',
                           ' -c ',
                           ' -bt ',
                           ' dodecahedron ',
                           ' -d ',
                           genbox_obj.get_box_size(),
                           ' -princ'])
        proc = Popen(command, shell=True, stdout=PIPE)
        proc.communicate()

        command = ''.join([gromacs_path.value,
                          './gmx',
                           ' solvate ',
                           ' -o ',
                           ' water ',
                           ' -cp ',
                           ' box.gro ',
                           ' -p ',
                           ' top.top ',
                           ' -cs ',
                           'spc216.gro'])
        proc = Popen(command, shell=True, stdout=PIPE)
        proc.communicate()

        command = ''.join([gromacs_path.value,
                          './gmx',
                           ' grompp ',
                           ' -f ',
                           min_none_mdp,
                           ' -c ',
                           ' water.gro ',
                           ' -p ',
                           ' top.top ',
                           ' -o ',
                           ' water.tpr'])
        proc = Popen(command, shell=True, stdout=PIPE)
        proc.communicate()

        command = ''.join([gromacs_path.value,
                           './gmx',
                           ' select ',
                           ' -f ',
                           ' water.gro ',
                           ' -s ',
                           ' water.tpr ',
                           ' -on ',
                           ' remover.ndx ',
                           ' -select ',
                           '\'"Agua remover" not same residue as resname SOL and (within 0.5 of group Protein)\''])
        proc = Popen(command, shell=True, stdout=PIPE)
        proc.communicate()

        command = ''.join([gromacs_path.value,
                           './gmx',
                           ' trjconv ',
                           ' -f ',
                           ' water.gro ',
                           ' -s ',
                           ' water.tpr ',
                           ' -n ',
                           ' remover.ndx ',
                           ' -o ',
                           ' water_ok.gro'])
        proc = Popen(command, shell=True, stdout=PIPE)
        proc.communicate()

        command = ''.join(['head -n-1',
                           ' top.top ',
                           '> top_temporary'])
        proc = Popen(command, shell=True, stdout=PIPE)
        proc.communicate()

        command = ''.join(['grep SOL',
                           ' water_ok.gro ',
                           '| grep OW | wc -l | awk \'{print $1}\''])
        proc = Popen(command, shell=True, stdout=PIPE)
        total_water, _ =  proc.communicate()
        total_water = total_water.decode('utf-8').rstrip()
        total_water = ''.join(['"',
                               total_water,
                               '"'])

        command = ''.join(['echo "SOL              "',
                          total_water,
                           '>> top_temporary'])
        proc = Popen(command, shell=True, stdout=PIPE)
        proc.communicate()

        command = ''.join(['mv top_temporary top.top'])
        proc = Popen(command, shell=True, stdout=PIPE)
        proc.communicate()

    list_obj_genbox = load_md_genbox(file_of_genbox)

    md_genboxRDD = sc.parallelize(list_obj_genbox)

    md_genboxRDD.foreach(run_genbox)

    finish_time = datetime.now()

    time_execution_log(finish_time, start_time, "gromacs_genbox.log")

