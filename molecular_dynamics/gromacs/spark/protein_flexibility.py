from pyspark import SparkContext, SparkConf, SparkFiles
from pyspark.sql import SQLContext, Row
import configparser
from subprocess import Popen, PIPE
from gromacs_utils import check_file_exists
import sys
import os
from os_utils import make_directory, preparing_path, time_execution_log
from pf_description import pf_description
from datetime import datetime


def load_prot_flex(file_of_protein_flexibility):
    list_ret = []
    f_file = open(file_of_protein_flexibility, "r")
    for line in f_file:
        splited_line = str(line).split()
        path = str(splited_line[0]).strip()
        nonwater_tpr = str(splited_line[1]).strip()
        echo_trjconv = str(splited_line[2]).strip()
        initial_frame_trjconv = int(splited_line[3])
        final_frame_trjconv = int(splited_line[4])
        nonwater_xtc = int(splited_line[5])
        echo_rmsf = str(splited_line[6]).strip()
        initial_frame_rmsf = int(splited_line[7])
        final_frame_rmsf = int(splited_line[8])
        title_output = str(splited_line[9]).strip()
        obj = pf_description(path,
                             nonwater_tpr,
                             echo_trjconv,
                             initial_frame_trjconv,
                             final_frame_trjconv,
                             nonwater_xtc,
                             echo_rmsf,
                             initial_frame_rmsf,
                             final_frame_rmsf,
                             title_output)
        list_ret.append(obj)

    return list_ret

def main():
    sc = SparkContext()
    sqlCtx = SQLContext(sc)

    config = configparser.ConfigParser()
    config.read('config.ini')

    # Path for gromacs spark project
    path_spark_gromacs = config.get('DRUGDESIGN', 'path_spark_gromacs')

    # Adding Python Source file
    sc.addPyFile(os.path.join(path_spark_gromacs, "gromacs_utils.py"))
    sc.addPyFile(os.path.join(path_spark_gromacs, "os_utils.py"))

    # Path for gromacs program
    gromacs_path = preparing_path(config.get('DRUGDESIGN', 'gromacs_path'))

    time_dt = int(config.get('GROMACS_ANALYSIS', 'time_dt'))

    # File that contains all md to create the trajectory
    file_of_prot_flex = sys.argv[1]
    check_file_exists(file_of_prot_flex)

    # Broadcast
    gromacs_path = sc.broadcast(gromacs_path)
    time_dt = sc.broadcast(time_dt)

    # Start time
    start_time = datetime.now()

# ********************* STARTING FUNCTION ***************************

    def run_protein_flexibility(pf_obj):
        flex_dir = os.path.join(pf_obj.get_path(), "flexibility")
        make_directory(flex_dir)

        nonwater_tpr = os.path.join(pf_obj.get_path(), pf_obj.get_nonwater_tpr())
        nonwater_xtc = os.path.join(pf_obj.get_path(), pf_obj.get_nonwater_xtc())

        traj_gro = os.path.join(flex_dir, "".join([pf_obj.get_title_output(),
                                                   ".",
                                                   pf_obj.get_output_sufix,
                                                   ".gro"]))
        rmsf_xvg = os.path.join(flex_dir, "".join(["rmsf_",
                                                   pf_obj.get_title_output(),
                                                   ".",
                                                   pf_obj.get_output_sufix,
                                                   ".xvg"]))


        command = "".join(["echo ",
                           pf_obj.get_echo_trjconv(),
                           " | ",
                           gromacs_path,
                           "./gmx trjconv",
                           " -s ",
                           nonwater_tpr,
                           " -f ",
                           nonwater_xtc,
                           " -o ",
                           traj_gro,
                           " -b ",
                           pf_obj.get_initial_frame_trjconv(),
                           " -e ",
                           pf_obj.get_final_frame_trjconv()])

        proc = Popen(command, shell=True, stdout=PIPE)
        proc.communicate()

        command = "".join(["echo ",
                           pf_obj.get_echo_rmsf(),
                           " | ",
                           gromacs_path,
                           "./gmx rmsf",
                           " -f ",
                           nonwater_xtc,
                           " -s ",
                           traj_gro,
                           " -o ",
                           rmsf_xvg,
                           " -b ",
                           pf_obj.get_initial_frame_rmsf(),
                           " -e ",
                           pf_obj.get_final_frame_rmsf(),
                           " -dt ",
                           time_dt,
                           " -res"])
        proc = Popen(command, shell=True, stdout=PIPE)
        proc.communicate()

    list_obj_pf = load_prot_flex(file_of_prot_flex)

    prot_flex = sc.parallelize(list_obj_pf)

    prot_flex.foreach(run_protein_flexibility)

    finish_time = datetime.now()

    time_execution_log(finish_time, start_time, "gromacs_flexibility.log")
