import ConfigParser as configparser
import os
import sys
from subprocess import Popen, PIPE
from pyspark import SparkContext, SparkConf, SparkFiles
from pyspark.sql import SQLContext, Row
from datetime import datetime
from os_utils import make_directory, preparing_path, time_execution_log, check_file_exists
from docking_description import vd_description


def load_vd_file(file_of_vina_docking):
    list_ret = []
    f_file = open(file_of_vina_docking, "r")
    next(f_file)
    for line in f_file:
        splited_line = str(line).split("\t")
        receptor = str(splited_line[0]).strip()
        ligand = str(splited_line[1]).strip()
        obj = vd_description(receptor,
                             ligand)
        list_ret.append(obj)

    return list_ret


if __name__ == '__main__':

    sc = SparkContext()
    sqlCtx = SQLContext(sc)

    config = configparser.ConfigParser()
    config.read('config.ini')

    # Vina configuration for broadcast
    config_vina = config.get('VINA', 'config_file')
    vina_path = config.get('VINA', 'vina_program')
    pdbqt_ligand_path = config.get('DEFAULT', 'pdbqt_ligand_path')
    pdbqt_receptor_path = config.get('DEFAULT', 'pdbqt_receptor_path')
    path_save_output = config.get('DEFAULT', 'path_save_structure')
    path_save_log = config.get('DEFAULT', 'path_save_log')
    path_spark_drugdesign = config.get('DRUGDESIGN', 'path_spark_drugdesign')

    path_save_log = preparing_path(path_save_log)
    make_directory(path_save_log)

    path_save_output = preparing_path(path_save_output)
    make_directory(path_save_output)

    # Adding Python Source file
    sc.addPyFile(os.path.join(path_spark_drugdesign, "docking_description.py"))

    # Broadcast
    vina_path = sc.broadcast(vina_path)
    pdbqt_ligand_path = sc.broadcast(pdbqt_ligand_path)
    pdbqt_receptor_path = sc.broadcast(pdbqt_receptor_path)
    path_save_output = sc.broadcast(path_save_output)
    path_save_log = sc.broadcast(path_save_log)
    sc.addFile(config_vina)

    file_of_vina_docking = sys.argv[1]
    check_file_exists(file_of_vina_docking)
    start_time = datetime.now()

    def run_vina_docking(vd_obj):

        receptor = ''.join([pdbqt_receptor_path.value,
                            vd_obj.get_receptor(),
                            ".pdbqt"])
        ligand = ''.join([pdbqt_ligand_path.value,
                          vd_obj.get_ligand(),
                          ".pdbqt"])
        output_save = ''.join([path_save_output.value,
                               vd_obj.get_receptor(),
                               "_-_",
                               vd_obj.get_ligand(),
                               ".pdbqt"])
        output_log = ''.join([path_save_log.value,
                              vd_obj.get_receptor(),
                              "_-_",
                              vd_obj.get_ligand(),
                              ".log"])

        command = "".join([vina_path.value,
                           " --config ",
                           config_vina,
                           " --receptor ",
                           receptor,
                           " --ligand ",
                           ligand,
                           " --out ",
                           output_save,
                           " --log ",
                           output_log])

        proc = Popen(command, shell=True, stdout=PIPE)
        proc.communicate()

    list_obj_vd = load_vd_file(file_of_vina_docking)
    vina_dockingRDD = sc.parallelize(list_obj_vd)
    vina_dockingRDD.foreach(run_vina_docking)

    finish_time = datetime.now()
    time_execution_log(finish_time, start_time, "vina_docking.log")