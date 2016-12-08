import configparser
import os
import sys
from pyspark import SparkContext, SparkConf, SparkFiles
from pyspark.sql import SQLContext, Row
from datetime import datetime
from subprocess import Popen, PIPE
from docking_description import sd_description
from os_utils import check_file_exists
from vina_utils import get_value_from_box_center, get_value_from_box_size, create_jsondata_from_docking_output_file, calculate_avg_value
import json


def load_sd_file(file_of_vina_docking):
    list_ret = []
    f_file = open(file_of_vina_docking, "r")
    for line in f_file:
        splited_line = str(line).split()
        pdb_id = str(splited_line[0]).strip()
        ligand = str(splited_line[1]).strip()
        obj = sd_description(pdb_id,
                           ligand)
        list_ret.append(obj)

    return list_ret

if __name__ == '__main__':

    sc = SparkContext()
    sqlCtx = SQLContext(sc)

    config = configparser.ConfigParser()
    config.read('config.ini')

    # Vina configuration for broadcast
    vina_program = config.get('VINA', 'vina_program')
    mol2_path = config.get('VINA', 'mol2_path')
    pythonsh = config.get('VINA', 'pythonsh')
    script_receptor4 = config.get('VINA', 'script_receptor4')
    script_ligand4 = config.get('VINA', 'script_ligand4')
    script_pdbqt_to_pdb = config.get('VINA', 'script_pdbqt_to_pdb')
    exhaustiveness = config.get('VINA', 'exhaustiveness')
    path_spark_drugdesign = config.get('DRUGDESIGN', 'path_spark_drugdesign')
    gromacs_path = config.get('DRUGDESIGN', 'gromacs_path')

    # Adding Python Source file
    sc.addPyFile(os.path.join(path_spark_drugdesign, "docking_description.py"))
    sc.addPyFile(os.path.join(path_spark_drugdesign, "os_utils.py"))
    sc.addPyFile(os.path.join(path_spark_drugdesign, "vina_utils.py"))


    # Broadcast
    vina_program = sc.broadcast(vina_program)
    pythonsh = sc.broadcast(pythonsh)
    script_receptor4 = sc.broadcast(script_receptor4)
    script_ligand4 = sc.broadcast(script_ligand4)
    script_pdbqt_to_pdb = sc.broadcast(script_pdbqt_to_pdb)
    exhaustiveness = sc.broadcast(exhaustiveness)
    gromacs_path = sc.broadcast(gromacs_path)

    file_of_single_docking = sys.argv[1]
    check_file_exists(file_of_single_docking)

    def run_single_docking(sing_dock_obj):

        ligand = ''.join([mol2_path, sing_dock_obj.get_ligand()])
        check_file_exists(ligand)
        docking_output_pdbqt = sing_dock_obj.get_pdb_id() + 'output.pdbqt'
        docking_output_log = sing_dock_obj.get_pdb_id() + '_output.log'


        # Downloads the receptor pdb file
        command = ''.join(["wget http://www.rcsb.org/pdb/files/",
                            sing_dock_obj.get_pdb_id(),
                           ".pdb >/dev/null 2>/dev/null"])
        proc = Popen(command, shell=True, stdout=PIPE)
        proc.communicate()

        # Removes waters, ligands and metals with ADT
        command = ''.join([pythonsh.value,
                           ' ',
                           script_receptor4.value,
                           ' -r ',
                           sing_dock_obj.get_pdb_id(),
                           '.pdb',
                           ' -o ',
                           ' temporary_',
                           sing_dock_obj.get_pdb_id(),
                          '.pdbqt ',
                           ' -A ',
                           ' none ',
                           ' -U ',
                           ' waters ',
                           ' -e ',
                           ' >/dev/null 2>/dev/null'])
        proc = Popen(command, shell=True, stdout=PIPE)
        proc.communicate()

        # Converts the .pdbqt back to .pdb
        command = ''.join([pythonsh.value,
                           ' ',
                           script_pdbqt_to_pdb.value,
                           ' -f ',
                           'temporary_',
                           sing_dock_obj.get_pdb_id(),
                           '.pdbqt',
                           ' -o ',
                           sing_dock_obj.get_pdb_id(),
                           '_ok.pdb',
                           '>/dev/null 2>/dev/null'])
        proc = Popen(command, shell=True, stdout=PIPE)
        proc.communicate()

        # Delete the temporary pdbqt
        command = ''.join(['rm temporary_',
                           sing_dock_obj.get_pdb_id(),
                          '.pdbqt'])
        proc = Popen(command, shell=True, stdout=PIPE)
        proc.communicate()

        # Uses gromacs to detect the box. note that the receptor is rotated!
        # This means that the pdb_id_box.pdb file must be used for visualizations

        command = ''.join([
                           'echo System | ',
                           gromacs_path.value,
                           'gmx editconf ',
                           ' -f ',
                           sing_dock_obj.get_pdb_id(),
                           '_ok.pdb',
                           ' -o ',
                           sing_dock_obj.get_pdb_id(),
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
                           sing_dock_obj.get_pdb_id(),
                           '_box.pdb'])
        proc = Popen(command, shell=True, stdout=PIPE)
        proc.communicate()

        command = ''.join([
                           "echo 'System' | ",
                           gromacs_path.value,
                           'gmx editconf ',
                           ' -f ',
                           sing_dock_obj.get_pdb_id(),
                           '_ok.pdb',
                           ' -o ',
                           sing_dock_obj.get_pdb_id(),
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
                           sing_dock_obj.get_pdb_id(),
                           '_ok.pdb'])
        proc = Popen(command, shell=True, stdout=PIPE)
        proc.communicate()

        # Prepare the ligand file
        command = ''.join([pythonsh.value,
                           ' ',
                           script_ligand4.value,
                           ' -l ',
                           ligand,
                           ' -o ',
                           'ligand.pdbqt ',
                           ' -A ',
                           ' none ',
                           ' -U ',
                           ' nphs_lps ',
                           '>/dev/null 2>/dev/null'])
        proc = Popen(command, shell=True, stdout=PIPE)
        proc.communicate()

        # Prepare the receptor file
        command = ''.join([pythonsh.value,
                           ' ',
                           script_receptor4.value,
                           ' -r ',
                           sing_dock_obj.get_pdb_id(),
                           '_box.pdb',
                           ' -o ',
                           'receptor.pdbqt ',
                           ' -A ',
                           ' none ',
                           ' -U ',
                           ' nphs_lps ',
                           '>/dev/null 2>/dev/null'])
        proc = Popen(command, shell=True, stdout=PIPE)
        proc.communicate()

        command = ''.join([vina_program.value,
                           ' --receptor ',
                           'receptor.pdbqt',
                           ' --ligand ',
                           'ligand.pdbqt',
                           ' --center_x ',
                           box_center.get('box_center_x', 'None'),
                           ' --center_y ',
                           box_center.get('box_center_y', 'None'),
                           ' --center_z ',
                           box_center.get('box_center_z', 'None'),
                           ' --size_x ',
                           box_size.get('box_size_x', 'None'),
                           ' --size_y ',
                           box_size.get('box_size_y', 'None'),
                           ' --size_z ',
                           box_size.get('box_size_z', 'None'),
                           ' --exhaustiveness ',
                           exhaustiveness.value,
                           ' --out ',
                           docking_output_pdbqt,
                           ' --log ',
                           docking_output_log,
                           '>/dev/null 2>/dev/null'])
        proc = Popen(command, shell=True, stdout=PIPE)
        proc.communicate()

        # Cleaning up
        command = ''.join(['rm ',
                           ' receptor.pdbqt ',
                           ' ligand.pdbqt'])
        proc = Popen(command, shell=True, stdout=PIPE)
        proc.communicate()

        # Creating Json output
        data = create_jsondata_from_docking_output_file(docking_output_log)
        output_file = sing_dock_obj.get_pdb_id() + '_log.json'
        with open(output_file, 'w') as outfile:
             json.dump(data, outfile, indent=4, sort_keys=True, separators=(',', ':'))
        data = calculate_avg_value(docking_output_log)
        output_file = sing_dock_obj.get_pdb_id() + '_avg.json'
        with open(output_file, 'w') as outfile:
             json.dump(data, outfile, indent=4, sort_keys=True, separators=(',', ':'))


    list_obj_single_dock = load_sd_file(file_of_single_docking)
    sing_dockRDD = sc.parallelize(list_obj_single_dock)
    sing_dockRDD.foreach(run_single_docking)
