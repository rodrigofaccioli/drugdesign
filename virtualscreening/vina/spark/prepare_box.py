import os
import sys
from pyspark import SparkContext, SparkConf, SparkFiles
from pyspark.sql import SQLContext, Row
from pbox_description import pbox_description
from subprocess import Popen, PIPE
from os_utils import check_file_exists, download_file
from vina_utils import get_value_from_box_center, get_value_from_box_size
from json_utils import create_json_file

try:
    import configparser
except ImportError:
    import ConfigParser as configparser


def load_pbox_file(file_of_pdbid_list):
    list_ret = []
    f_file = open(file_of_pdbid_list, "r")
    for line in f_file:
        splited_line = str(line).split()
        pdb_id = str(splited_line[0]).strip()
        obj = pbox_description(pdb_id)
        list_ret.append(obj)

    return list_ret

if __name__ == '__main__':

    sc = SparkContext()
    sqlCtx = SQLContext(sc)

    config = configparser.ConfigParser()
    config.read('config.ini')

    # Vina configuration for broadcast
    vina_program = config.get('VINA', 'vina_program')    
    pythonsh = config.get('VINA', 'pythonsh')
    script_receptor4 = config.get('VINA', 'script_receptor4')
    script_ligand4 = config.get('VINA', 'script_ligand4')
    script_pdbqt_to_pdb = config.get('VINA', 'script_pdbqt_to_pdb')
    exhaustiveness = config.get('VINA', 'exhaustiveness')
    path_spark_drugdesign = config.get('DRUGDESIGN', 'path_spark_drugdesign')
    gromacs_path = config.get('DRUGDESIGN', 'gromacs_path')

    # Adding Python Source file
    sc.addPyFile(os.path.join(path_spark_drugdesign, "pbox_description.py"))
    sc.addPyFile(os.path.join(path_spark_drugdesign, "os_utils.py"))
    sc.addPyFile(os.path.join(path_spark_drugdesign, "vina_utils.py"))
    sc.addPyFile(os.path.join(path_spark_drugdesign, "json_utils.py"))


    # Broadcast
    vina_program = sc.broadcast(vina_program)
    pythonsh = sc.broadcast(pythonsh)
    script_receptor4 = sc.broadcast(script_receptor4)
    script_ligand4 = sc.broadcast(script_ligand4)
    script_pdbqt_to_pdb = sc.broadcast(script_pdbqt_to_pdb)
    exhaustiveness = sc.broadcast(exhaustiveness)
    gromacs_path = sc.broadcast(gromacs_path)

    file_of_pdbid_list = sys.argv[1]
    check_file_exists(file_of_pdbid_list)

    def run_prepare_box(prepare_box_obj):


        # Downloads the receptor pdb file
        file_name = prepare_box_obj.get_pdb_id() + '.pdb'
        url = 'http://www.rcsb.org/pdb/files/' + file_name
        download_file(url, file_name)

        # Removes waters, ligands and metals with ADT
        command = ''.join([pythonsh.value,
                           ' ',
                           script_receptor4.value,
                           ' -r ',
                           prepare_box_obj.get_pdb_id(),
                           '.pdb',
                           ' -o ',
                           'temporary_',
                           prepare_box_obj.get_pdb_id(),
                           '.pdbqt ',
                           '-A ',
                           'none ',
                           '-U ',
                           'waters ',
                           '-e ',
                           '>/dev/null 2>/dev/null'])
        proc = Popen(command, shell=True, stdout=PIPE)
        proc.communicate()

        # Converts the .pdbqt back to .pdb
        command = ''.join([pythonsh.value,
                           ' ',
                           script_pdbqt_to_pdb.value,
                           ' -f ',
                           'temporary_',
                           prepare_box_obj.get_pdb_id(),
                           '.pdbqt',
                           ' -o ',
                           prepare_box_obj.get_pdb_id(),
                           '_ok.pdb',
                           '>/dev/null 2>/dev/null'])
        proc = Popen(command, shell=True, stdout=PIPE)
        proc.communicate()


        # Delete the temporary pdbqt
        command = ''.join(['rm temporary_',
                           prepare_box_obj.get_pdb_id(),
                           '.pdbqt'])
        proc = Popen(command, shell=True, stdout=PIPE)
        proc.communicate()

        # Uses gromacs to detect the box. note that the receptor is rotated!
        # This means that the pdb_id_box.pdb file must be used for visualizations

        command = ''.join(['echo System | ',
                           gromacs_path.value,
                           'gmx editconf ',
                           ' -f ',
                           prepare_box_obj.get_pdb_id(),
                           '_ok.pdb',
                           ' -o ',
                           prepare_box_obj.get_pdb_id(),
                           '_box.pdb',
                           ' -d ',
                           ' 0.5 ',
                           ' -bt ',
                           ' triclinic ',
                           ' -angles ',
                           ' 90 90 90 ',
                           ' -princ ',
                           ' -c ',
                           ' 2>/dev/null '
                           "| grep 'new center' ",
                           "| awk '{print $4,$5,$6}' "])
        proc = Popen(command, shell=True, stdout=PIPE)
        box_center, _ = proc.communicate()
        box_center = box_center.decode("utf-8")
        box_center = get_value_from_box_center(box_center)

        command = ''.join(['rm ',
                           prepare_box_obj.get_pdb_id(),
                           '_box.pdb'])
        proc = Popen(command, shell=True, stdout=PIPE)
        proc.communicate()

        command = ''.join(["echo 'System' | ",
                           gromacs_path.value,
                           'gmx editconf ',
                           ' -f ',
                           prepare_box_obj.get_pdb_id(),
                           '_ok.pdb',
                           ' -o ',
                           prepare_box_obj.get_pdb_id(),
                           '_box.pdb',
                           ' -d ',
                           ' 0.5 ',
                           ' -bt ',
                           ' triclinic ',
                           ' -angles ',
                           ' 90 90 90 ',
                           ' -princ ',
                           ' -c ',
                           '2>/dev/null',
                           "| grep 'new box vectors' ",
                           "| awk '{print $5,$6,$7}' "])
        proc = Popen(command, shell=True, stdout=PIPE)
        box_size, _ = proc.communicate()
        box_size = box_size.decode("utf-8")
        box_size = get_value_from_box_size(box_size)

        # Deletes the temporary .pdb
        command = ''.join(['rm ',
                           prepare_box_obj.get_pdb_id(),
                           '.pdb'])
        proc = Popen(command, shell=True, stdout=PIPE)
        proc.communicate()

        command = ''.join(['rm ',
                           prepare_box_obj.get_pdb_id(),
                           '_box.pdb'])
        proc = Popen(command, shell=True, stdout=PIPE)
        proc.communicate()

        command = ''.join(['rm ',
                           prepare_box_obj.get_pdb_id(),
                           '_ok.pdb'])
        proc = Popen(command, shell=True, stdout=PIPE)
        proc.communicate()

        # Creating Json output

        data = dict(box_center, **box_size)
        output_file_name = prepare_box_obj.get_pdb_id() + '_box.json'
        create_json_file(output_file_name, data)


    list_obj_prepare_box = load_pbox_file(file_of_pdbid_list)
    sing_dockRDD = sc.parallelize(list_obj_prepare_box)
    sing_dockRDD.foreach(run_prepare_box)
