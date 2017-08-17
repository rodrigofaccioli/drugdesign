import os
import ntpath
import json
from math import sqrt
#from json_utils import create_json_file


def get_files_mol2(mypath):
    only_mol2_file = []
    for root, dirs, files in os.walk(mypath):
        for file in files:
            if file.endswith(".mol2"):
                f_path = os.path.join(root, file)
                only_mol2_file.append(f_path)
    return only_mol2_file


""" This function obtains all pdb files
    in mypath
"""


def get_files_pdb(mypath):
    only_pdb_file = []
    for root, dirs, files in os.walk(mypath):
        for file in files:
            if file.endswith(".pdb"):
                f_path = os.path.join(root, file)
                only_pdb_file.append(f_path)
    return only_pdb_file


""" This function obtains all pdb files
    in mypath filtered by reference
"""


def get_files_pdb_filter(mypath, reference):
    only_pdb_file = []
    for root, dirs, files in os.walk(mypath):
        for file in files:
            if file.endswith(".pdb"):
                if file.find(reference) > -1:
                    f_path = os.path.join(root, file)
                    only_pdb_file.append(f_path)
    return only_pdb_file


""" This function obtains all pdbqt files
    in mypath
"""


def get_files_pdbqt(mypath):
    only_pdb_file = []
    for root, dirs, files in os.walk(mypath):
        for file in files:
            if file.endswith(".pdbqt"):
                f_path = os.path.join(root, file)
                only_pdb_file.append(f_path)
    return only_pdb_file


""" This function obtains all log files
    in mypath
"""


def get_files_log(mypath):
    only_pdb_file = []
    for root, dirs, files in os.walk(mypath):
        for file in files:
            if file.endswith(".log"):
                f_path = os.path.join(root, file)
                only_pdb_file.append(f_path)
    return only_pdb_file


""" This function obtains the name of
sorted energy file
"""


def get_file_name_sorted_energy():
    return 'vs_energies_sorted.txt'


def get_separator_filename_mode():
    """
	Returns the separator file mode. It means a way to separate receptor_ligand and mode
    Example:
        >>> get_separator_filename_mode()
    @return: the separator file mode
    @rtype: string
	"""
    return '+----+'


""" This function obtains the name of
path that saving pdbqt files for analysis
"""


def get_directory_pdbqt_analysis(path_analysis):
    path_analysis_pdbqt = os.path.join(path_analysis, "pdbqt_model")
    # Checking path_analysis
    if not os.path.exists(path_analysis_pdbqt):
        os.makedirs(path_analysis_pdbqt)
    return path_analysis_pdbqt


def get_structure_file_name(myfile):
    """
	This function obtains the name of myfile without filename extension
	"""
    path, filename = ntpath.split(myfile)
    name = str(filename.split(".")[0])  # remove .pdbqt
    return name


def get_name_model_pdb(myfile):
    """
	This function obtains the name of myfile without filename extension
	"""
    path, filename = ntpath.split(myfile)
    name = str(filename.split(".")[0])  # remove .pdb
    return name


def get_name_model_pdbqt(myfile):
    """
	This function obtains the name of myfile without filename extension
	"""
    path, filename = ntpath.split(myfile)
    name = str(filename.split(".")[0])  # remove .pdbqt
    return name


def get_name_receptor_pdb(myfile):
    """
	This function obtains the name of myfile without filename extension
	"""
    path, filename = ntpath.split(myfile)
    name = str(filename.split(".")[0])  # remove .pdb
    return name


def get_name_receptor_pdbqt(myfile):
    """
	This function obtains the name of myfile without filename extension
	"""
    path, filename = ntpath.split(myfile)
    name = str(filename.split(".")[0])  # remove .pdbqt
    return name


def get_ligand_from_receptor_ligand_model(receptor_ligand_model):
    """
	This function obtains the name of ligand based on receptor_ligand_model
	Example of input: compl_ns3pro_dm_0_-_NuBBE_485_obabel_3D+----+20
	"""
    separator_model = get_separator_filename_mode()
    separator_receptor = "_-_"
    string_ref = receptor_ligand_model

    s = string_ref.split(separator_receptor)  # Removing receptor
    s = str(s[1]).split(separator_model)  # Removing model
    ligand_name = str(s[0])  # geting name of ligand
    return ligand_name


""" This function obtains the name of
path that saving pdbqt files for analysis
"""


def get_directory_pdb_analysis(path_analysis):
    path_analysis_pdb = os.path.join(path_analysis, "pdb_model")
    # Checking path_analysis
    if not os.path.exists(path_analysis_pdb):
        os.makedirs(path_analysis_pdb)
    return path_analysis_pdb


""" This function obtains the name of
path that saving pdbqt files for analysis
"""


def get_directory_complex_pdb_analysis(path_analysis):
    path_analysis_pdb = os.path.join(path_analysis, "pdb_complex")
    # Checking path_analysis
    if not os.path.exists(path_analysis_pdb):
        os.makedirs(path_analysis_pdb)
    return path_analysis_pdb


""" This function loading pdb file to list.
list_ret is composed by pdb_path_file and loaded file.
"""


def loading_pdb_2_list(pdb_path_file):
    list_pdb = []
    f_PDB = open(pdb_path_file, "r")
    for line in f_PDB:
        if line.find("ATOM") > -1:
            list_pdb.append(line)
    f_PDB.close()
    list_ret = (pdb_path_file, list_pdb)
    return list_ret


"""	Extract the numbers. GROMACS uses nm and autodock vina uses A as units.
Therefore, values from gromacs are multiplied by 10.
"""


def get_value_from_box_center(box):
    splited_value_box = str(box).split()
    return dict(box_center_x="{0:.2f}".format(float(splited_value_box[0]) * 10),
                box_center_y="{0:.2f}".format(float(splited_value_box[1]) * 10),
                box_center_z="{0:.2f}".format(float(splited_value_box[2]) * 10))


def get_value_from_box_size(box):
    splited_value_box = str(box).split()
    return dict(box_size_x="{0:.2f}".format(float(splited_value_box[0]) * 10),
                box_size_y="{0:.2f}".format(float(splited_value_box[1]) * 10),
                box_size_z="{0:.2f}".format(float(splited_value_box[2]) * 10))


def _generate_parameters_to_complexo_dm():
    """These parameters are generated in hardcoded form, so this function will be deprecated later"""
    d = dict(num_modes=9999,
             energy_range=9999,
             exhaustiveness=10,
             cpu=1)
    create_json_file('parameters_complexo_dm.json', d)


# Use pdbid_box.json and general_parameters.json to crete the config_complexo_dm.txt
def generate_config_complexo_dm(box_json, general_parameters_json):
    box = open(box_json, 'r')
    box_dict = json.load(box)
    parameters = open(general_parameters_json, 'r')
    parameters_dict = json.load(parameters)

    file = open('config_complexo_dm.txt', 'w')

    for key, value in box_dict.items():
        line = ''.join([str(key),
                        ' = ',
                        str(value),
                        '\n'])
        file.write(line)

    for key, value in parameters_dict.items():
        line = ''.join([str(key),
                        ' = ',
                        str(value),
                        '\n'])
        file.write(line)
    file.close()


def calculate_avg_value(docking_output):
    f_file = open(docking_output, "r")
    err_list = []
    value = 0
    n = 0
    for line in f_file:
        splited_line = str(line).split()
        if not len(splited_line) == 0 and splited_line[0].isdigit():
            value += float(splited_line[1])
            err_list.append(float(splited_line[2]))
            n += 1

    avg = value / n
    avg = "%.1f" % avg
    err = 0

    for value in err_list:
        err = (value - float(avg)) * (value - float(avg))
        err = sqrt(err / (n - 1) / sqrt(n))
        err = "%.1f" % err

    f_file.close()
    return dict(number_modes=n,
                avg=avg,
                err=err)

def get_receptor_from_receptor_ligand_model(receptor_ligand_model):
    """
    This function obtains the name of receptor based on receptor_ligand_model
    Example of input: compl_ns3pro_dm_0_-_NuBBE_485_obabel_3D+----+20
    """
    separator_model = get_separator_filename_mode()
    separator_receptor = "_-_"
    string_ref = receptor_ligand_model

    receptor_name = string_ref.split(separator_receptor)[0] #Removing all, except receptor name
    return receptor_name

def get_model_from_receptor_ligand_model(receptor_ligand_model):
    """
    This function obtains the model based on receptor_ligand_model
    Example of input: compl_ns3pro_dm_0_-_NuBBE_485_obabel_3D+----+20
    Return: 20
    """
    separator_model = get_separator_filename_mode()
    separator_receptor = "_-_"
    string_ref = receptor_ligand_model

    s = string_ref.split(separator_receptor) #Removing receptor
    s = str(s[1]).split(separator_model) #Removing ligand name
    # getting model
    model = int(s[len(s)-1])
    return model

""" This function obtains the name of
path that saving pdbqt files for analysis
"""
def get_directory_temp_analysis(path_analysis):
    path_analysis_temp = os.path.join(path_analysis,"temp")
    #Checking path_analysis
    if not os.path.exists(path_analysis_temp):
        os.makedirs(path_analysis_temp)
    return path_analysis_temp
