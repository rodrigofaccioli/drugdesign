import os
from vina_utils import get_separator_filename_mode, get_structure_file_name, get_name_model_pdb
from subprocess import Popen, PIPE

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



def pdbqt2pdb(param):
	model = str(param[0])
	path_analysis_pdb = str(param[1])
	pythonsh = str(param[2])
	pdbqt_to_pdb = str(param[3])

	fpdb_filename = get_name_model_pdb(model)
	fpdb_filename = fpdb_filename+".pdb"
	fpdb = os.path.join(path_analysis_pdb,fpdb_filename)				
	process = Popen([pythonsh, pdbqt_to_pdb, '-f', model, '-o', fpdb], stdout=PIPE, stderr=PIPE)		
	stdout, stderr = process.communicate()

""" This function obtains all atoms that are
in input list. It returns a list of all atom
"""
def get_atom_section_from_atom_list(path_file_name, atom_list):
	return_list = []
	f_file = open(path_file_name, "r")
	for line in f_file:
		if (line.find("ATOM") > -1): #| (line.find("HETATM") > -1)
			return_list.append(line.rstrip())
#			splited_line = line.split()
#			last_collum_from_line = str(splited_line[-1])#Get the last column
#			for atom in atom_list:
#				if last_collum_from_line.find(atom) > -1: #found atom
#					return_list.append(line.rstrip())				 
	f_file.close()
	return return_list

""" This function save a text file from List
"""
def save_text_file_from_list(path_file_for_saving, list_of_tupla_ref):
	f_file = open(path_file_for_saving,"w")
	for item in list_of_tupla_ref:
		line = ""
		for i in item:
			line += str(i) +" "
		line += "\n" 
		f_file.write(line)
	f_file.close()
