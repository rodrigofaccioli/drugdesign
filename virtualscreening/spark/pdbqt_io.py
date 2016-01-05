import os
from vina_utils import get_separator_filename_mode, get_structure_file_name

def save_pdbqt_from_list(list_line, path_save, base_file_name_model):
	pdbqt_file = os.path.join(path_save, base_file_name_model)
	f_file = open(pdbqt_file, "w")
	for item in list_line:
		f_file.write(item)
	f_file.close()


def split_pdbqt(param):
	list_line = []
	structure = str(param[0])
	path_save = str(param[1])

	#Preparing base file name
	#Obtaing base file name from structure
	base_file_name = get_structure_file_name(structure)
	#Adding model separator in base_file_name
	base_file_name = base_file_name + get_separator_filename_mode()

	#Loading and spliting models in pdbqt files	
	f_file = open(structure, "r")
	for line in f_file:
		list_line.append(line)
		if line.find("MODEL") > -1:
			current_model = line.split()[1]
			base_file_name_model = base_file_name+current_model+".pdbqt"
		elif line.find("ENDMDL") > -1: #Saving current model
			save_pdbqt_from_list(list_line, path_save, base_file_name_model)
			list_line = []		



