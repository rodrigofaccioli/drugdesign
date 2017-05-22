from pyspark import SparkContext, SparkConf, SparkFiles
from pyspark.sql import SQLContext, Row
from subprocess import Popen, PIPE
from datetime import datetime
import os
import sys
import configparser
from md_description import md_description
from gromacs_utils import check_file_exists
from os_utils import make_directory, preparing_path, time_execution_log
from basic_analysis import run_basic_analysis


def search_for_ligand_ndx_file(ndx_file):
    other = 0
    textFile = open(ndx_file, "r")
    for line in textFile:
        if line.upper().find("OTHER") > -1:
            other = other + 1
    textFile.close()
    # Count how many other there is
    if other > 0:
        return True
    else:
        return False


def load_md_traj(file_of_md_analysis):
    list_ret = []
    f_file = open(file_of_md_analysis, "r")
    for line in f_file:
        splited_line = str(line).split()
        path = str(splited_line[0]).strip()
        prefix_ref = str(splited_line[1]).strip()
        repetion_number = int(splited_line[2])
        title_output = str(splited_line[3]).strip()
        total_running = int(splited_line[4])
        obj = md_description(path,
                             prefix_ref,
                             repetion_number,
                             title_output,
                             total_running)
        list_ret.append(obj)

    return list_ret


def main():
    sc = SparkContext()
    sqlCtx = SQLContext(sc)

    config = configparser.ConfigParser()
    config.read('config.ini')

    # Path for gromacs spark project
    path_spark_drugdesign = config.get('DRUGDESIGN', 'path_spark_drugdesign')

    # Adding Python Source file
    sc.addPyFile(os.path.join(path_spark_drugdesign, "gromacs_utils.py"))
    sc.addPyFile(os.path.join(path_spark_drugdesign, "os_utils.py"))
    sc.addPyFile(os.path.join(path_spark_drugdesign, "basic_analysis.py"))
    sc.addPyFile(os.path.join(path_spark_drugdesign, "md_description.py"))

    # Path for gromacs program
    gromacs_path = preparing_path(config.get('DRUGDESIGN', 'gromacs_path'))

    time_dt = int(config.get('GROMACS_ANALYSIS', 'time_dt'))
    time_dt_pdb = int(config.get('GROMACS_ANALYSIS', 'time_dt_pdb'))
    water_layer_thickness = int(config.get('GROMACS_ANALYSIS',
                                           'water_layer_thickness'))

    # File that contains all md to create the trajectory
    file_of_md_analysis = sys.argv[1]
    check_file_exists(file_of_md_analysis)

    start_time = datetime.now()

    # Broadcast
    gromacs_path = sc.broadcast(gromacs_path)
    time_dt = sc.broadcast(time_dt)
    time_dt_pdb = sc.broadcast(time_dt_pdb)
    water_layer_thickness = sc.broadcast(water_layer_thickness)

# ********************* STARTING FUNCTION ***************************
    def run_trajetory(md_obj):
        ana_dir = os.path.join(md_obj.get_path(), "analysis")
        make_directory(ana_dir)

        # Original file names from the simulation
        reference_xtc = os.path.join(md_obj.get_path(),
                                     md_obj.get_simulation_prefix() + ".xtc")
        reference_tpr = os.path.join(md_obj.get_path(),
                                     md_obj.get_simulation_prefix() + ".tpr")

        # File names after trajectory treatment.
        allatom_xtc = os.path.join(ana_dir, "".join([md_obj.get_prefix_ref(),
                                                     "_fit.",
                                                     str(md_obj.get_repetion_number()),
                                                     ".xtc"]))
        allatom_tpr = reference_tpr
        nonwater_xtc = os.path.join(ana_dir,"".join([md_obj.get_prefix_ref(),
                                                     "_non-water.",
                                                     str(md_obj.get_repetion_number()),
                                                     ".xtc"]))
        nonwater_tpr = os.path.join(ana_dir, "".join([md_obj.get_prefix_ref(),
                                                      "_non-water.",
                                                      str(md_obj.get_repetion_number()),
                                                      ".tpr"]))
        nonwater_pdb = os.path.join(ana_dir, "".join([md_obj.get_prefix_ref(),
                                                      "_non-water.",
                                                      str(md_obj.get_repetion_number()),
                                                      ".pdb"]))
        waterlayer_pdb = os.path.join(ana_dir, "".join([md_obj.get_prefix_ref(),
                                                        "_water-",
                                                        str(water_layer_thickness.value),
                                                        "A-layer.",
                                                        str(md_obj.get_repetion_number()),
                                                        ".pdb"]))

        # Trajectory treatment to remove PBC artifacts
        xtc_whole = os.path.join(ana_dir, "".join([md_obj.get_prefix_ref(),
                                                   "_whole.",
                                                   str(md_obj.get_repetion_number()),
                                                   ".xtc"]))

        command = "".join(["echo System | ",
                           gromacs_path.value,
                           "./gmx trjconv ",
                           "-f ",
                           reference_xtc,
                           " -s ",
                           reference_tpr,
                           " -pbc whole",
                           " -o ",
                           xtc_whole,
                           " >/dev/null 2>/dev/null"])
        proc = Popen(command, shell=True, stdout=PIPE)
        proc.communicate()

        # Extracting first frame
        gro_first_frame = os.path.join(ana_dir, "".join(["0.",
                                                         str(md_obj.get_repetion_number()),
                                                         ".gro"]))
        command = "".join(["echo System | ",
                           gromacs_path.value,
                           "./gmx trjconv ",
                           "-f ",
                           xtc_whole,
                           " -s ",
                           reference_tpr,
                           " -e 0.1 ",
                           " -o ",
                           gro_first_frame,
                           " >/dev/null 2>/dev/null"])
        proc = Popen(command, shell=True, stdout=PIPE)
        proc.communicate()

        # Removing jumps
        xtc_nojump = os.path.join(ana_dir,
                                  "".join([md_obj.get_prefix_ref(),
                                           "_nojump.",
                                           str(md_obj.get_repetion_number()),
                                           ".xtc"]))
        command = "".join(["echo System | ",
                           gromacs_path.value,
                           "./gmx trjconv ",
                           "-f ",
                           xtc_whole,
                           " -s ",
                           gro_first_frame,
                           " -pbc nojump ",
                           " -o ",
                           xtc_nojump,
                           " >/dev/null 2>/dev/null"])
        proc = Popen(command, shell=True, stdout=PIPE)
        proc.communicate()

        # Centering the protein
        xtc_center_protein = os.path.join(ana_dir, "".join([md_obj.get_prefix_ref(),
                                                            "_center.",
                                                            str(md_obj.get_repetion_number()),
                                                            ".xtc"]))
        command = "".join(["echo C-alpha System | ",
                           gromacs_path.value,
                           "./gmx trjconv ",
                           "-f ",
                           xtc_whole,
                           " -s ",
                           gro_first_frame,
                           " -center ",
                           " -o ",
                           xtc_center_protein,
                           " >/dev/null 2>/dev/null"])
        proc = Popen(command, shell=True, stdout=PIPE)
        proc.communicate()

        # Putting all atoms in a compact box
        xtc_atoms_box = os.path.join(ana_dir, "".join([md_obj.get_prefix_ref(),
                                                       "_atom.",
                                                       str(md_obj.get_repetion_number()),
                                                       ".xtc"]))
        command = "".join(["echo System | ",
                           gromacs_path.value,
                           "./gmx trjconv ",
                           "-f ",
                           xtc_center_protein,
                           " -s ",
                           gro_first_frame,
                           " -ur compact ",
                           " -pbc atom ",
                           " -o ",
                           xtc_atoms_box,
                           " >/dev/null 2>/dev/null"])
        proc = Popen(command, shell=True, stdout=PIPE)
        proc.communicate()

        # Fitting the protein
        command = "".join(["echo C-alpha System | ",
                           gromacs_path.value,
                           "./gmx trjconv ",
                           "-f ",
                           xtc_atoms_box,
                           " -s ",
                           gro_first_frame,
                           " -fit rot+trans ",
                           " -o ",
                           allatom_xtc,
                           " >/dev/null 2>/dev/null"])
        proc = Popen(command, shell=True, stdout=PIPE)
        proc.communicate()

        # Creating water-free trajectory
        command = "".join(["echo non-water | ",
                           gromacs_path.value,
                           "./gmx convert-tpr ",
                           " -s ",
                           reference_tpr,
                           " -o ",
                           nonwater_tpr,
                           " >/dev/null 2>/dev/null"])
        proc = Popen(command, shell=True, stdout=PIPE)
        proc.communicate()
        command = "".join(["echo non-water | ",
                           gromacs_path.value,
                           "./gmx trjconv ",
                           "-f ",
                           allatom_xtc,
                           " -s ",
                           gro_first_frame,
                           " -o ",
                           nonwater_xtc,
                           " >/dev/null 2>/dev/null"])
        proc = Popen(command, shell=True, stdout=PIPE)
        proc.communicate()
        command = "".join(["echo system | ",
                           gromacs_path.value,
                           "./gmx trjconv ",
                           " -f ",
                           nonwater_xtc,
                           " -s ",
                           nonwater_tpr,
                           " -o ",
                           nonwater_pdb,
                           " -dt ",
                           str(time_dt_pdb.value),
                           " >/dev/null 2>/dev/null"])
        proc = Popen(command, shell=True, stdout=PIPE)
        proc.communicate()

        # Creating water_layer_thickness - A water-layer pdb trajectory
        t = 0
        frame = 0
        ndx_water_layer = os.path.join(ana_dir, "".join([md_obj.get_prefix_ref(),
                                                         "_water-layer.",
                                                         str(md_obj.get_repetion_number()),
                                                         ".ndx"]))
        ndx_temporary = os.path.join(ana_dir, "".join([md_obj.get_prefix_ref(),
                                                       "_temporary_",
                                                       str(md_obj.get_repetion_number()),
                                                       ".ndx"]))
        if os.path.isfile(waterlayer_pdb):
            os.remove(waterlayer_pdb)
        if os.path.isfile(ndx_water_layer):
            os.remove(ndx_water_layer)
        select_string = ('\'"water_layer" (same residue as ((resname SOL and within 0.'"$water_layer_thickness"' of group "Protein"))) or\
                        (group "Ion" and within 0.'"$water_layer_thickness"' of group "Protein") \
                         or (group "Protein") \'')
        select_string = select_string.replace("$water_layer_thickness",
                                              str(water_layer_thickness.value))
        # Running make_ndx
        command = "".join(["echo -e ",
                           "'chain z'\"\\n\"'q'\"\\n\" | ",
                           gromacs_path.value,
                           "gmx make_ndx ",
                           "-f ",
                           reference_tpr,
                           " -o ",
                           ndx_temporary,
                           " >/dev/null 2>/dev/null"])
        proc = Popen(command, shell=True, stdout=PIPE)
        proc.communicate()

        # Are there ligands?
        if search_for_ligand_ndx_file(ndx_temporary) is True:
            select_string = (select_string
                             + '\'or (same residue as ((resname SOL and within 0.'"$water_layer_thickness"' of group "Other"))) \
                             or (group "Ion" and within 0.'"$water_layer_thickness"' of group "Other") \
                             or (group "Other")\'')
        select_string = select_string.replace("$water_layer_thickness",
                                              str(water_layer_thickness.value))
        command = "".join([gromacs_path.value,
                           "gmx select -f ",
                           allatom_xtc,
                           " -s ",
                           allatom_tpr,
                           " -on ",
                           ndx_water_layer,
                           " -select ",
                           select_string,
                           " -dt ",
                           str(time_dt_pdb.value),
                           " >/dev/null 2>/dev/null"])
        proc = Popen(command, shell=True, stdout=PIPE)
        proc.communicate()

        # Creating pdb files
        command = "".join(["echo ",
                           str(frame),
                           " | ",
                           gromacs_path.value,
                           "./gmx trjconv ",
                           "-f ",
                           allatom_xtc,
                           " -s ",
                           allatom_tpr,
                           " -n ",
                           ndx_water_layer,
                           " -o ",
                           "frame_",
                           str(frame),
                           ".pdb ",
                           "-b ",
                           str(t),
                           " -e ",
                           str(t),
                           " >/dev/null 2>/dev/null"])
        proc = Popen(command, shell=True, stdout=PIPE)
        proc.communicate()
        command = "".join(["echo MODEL ", str(frame), " >> ", waterlayer_pdb])
        proc = Popen(command, shell=True, stdout=PIPE)
        proc.communicate()
        command = "".join(["grep ATOM ",
                           "frame_",
                           str(frame),
                           ".pdb ",
                           ">> ",
                           waterlayer_pdb])
        proc = Popen(command, shell=True, stdout=PIPE)
        proc.communicate()
        command = "".join(["echo ENDML", ">> ", waterlayer_pdb])
        proc = Popen(command, shell=True, stdout=PIPE)
        proc.communicate()

        # Removing temporary files
        command = "".join(["rm frame_", str(frame), ".pdb"])
        proc = Popen(command, shell=True, stdout=PIPE)
        proc.communicate()
        frame = frame + 1
        t = t + int(time_dt_pdb.value)

        if os.path.isfile(xtc_whole):
            os.remove(xtc_whole)
        if os.path.isfile(xtc_nojump):
            os.remove(xtc_nojump)
        if os.path.isfile(xtc_center_protein):
            os.remove(xtc_center_protein)
        if os.path.isfile(xtc_atoms_box):
            os.remove(xtc_atoms_box)
        if os.path.isfile(ndx_water_layer):
            os.remove(ndx_water_layer)
        if os.path.isfile(gro_first_frame):
            os.remove(gro_first_frame)
        command = "rm \#* 2>/dev/null"
        proc = Popen(command, shell=True, stdout=PIPE)
        proc.communicate()

        # Basic Analysis
        basic_an_data = (gromacs_path.value,
                         nonwater_xtc,
                         nonwater_tpr,
                         md_obj.get_simulation_prefix(),
                         ana_dir,
                         time_dt.value)
        run_basic_analysis(basic_an_data)

# ************************** END FUNCTION **********************************

    list_obj_md = load_md_traj(file_of_md_analysis)

    md_trajRDD = sc.parallelize(list_obj_md)

    md_trajRDD.foreach(run_trajetory)

    finish_time = datetime.now()

    time_execution_log(finish_time, start_time, "gromacs_trajectory.log")

main()
