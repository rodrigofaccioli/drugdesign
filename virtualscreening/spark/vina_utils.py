import os
import ntpath

def get_files_mol2(mypath):
	only_mol2_file = []
	for root, dirs, files in os.walk(mypath):
		for file in files:
			if file.endswith(".mol2"):
				f_path = os.path.join(root,file)
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
				f_path = os.path.join(root,file)
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
				f_path = os.path.join(root,file)
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
				f_path = os.path.join(root,file)
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
	path_analysis_pdbqt = os.path.join(path_analysis,"pdbqt_model")
	#Checking path_analysis
	if not os.path.exists(path_analysis_pdbqt):
		os.makedirs(path_analysis_pdbqt)
	else:
		if len(os.listdir(path_analysis_pdbqt)) > 0:
			raise EnvironmentError("Analysis directory for PDBQT contains files ")
	return path_analysis_pdbqt

def get_structure_file_name(myfile):
	""" 
	This function obtains the name of myfile without filename extension
	"""	
	path, filename = ntpath.split(myfile)
	name =  str(filename.split(".")[0]) #remove .pdbqt
	return name

def get_name_model_pdb(myfile):
	""" 
	This function obtains the name of myfile without filename extension
	"""	
	path, filename = ntpath.split(myfile)
	name =  str(filename.split(".")[0]) #remove .pdb
	return name

def get_name_receptor_pdb(myfile):
	""" 
	This function obtains the name of myfile without filename extension
	"""	
	path, filename = ntpath.split(myfile)
	name =  str(filename.split(".")[0]) #remove .pdb
	return name

""" This function obtains the name of 
path that saving pdbqt files for analysis 
"""
def get_directory_pdb_analysis(path_analysis):
	path_analysis_pdb = os.path.join(path_analysis,"pdb_model")
	#Checking path_analysis
	if not os.path.exists(path_analysis_pdb):
		os.makedirs(path_analysis_pdb)
	else:
		if len(os.listdir(path_analysis_pdb)) > 0:
			raise EnvironmentError("Analysis directory for PDB contains files ")
	return path_analysis_pdb

""" This function obtains the name of 
path that saving pdbqt files for analysis 
"""
def get_directory_complex_pdb_analysis(path_analysis):
	path_analysis_pdb = os.path.join(path_analysis,"pdb_complex")
	#Checking path_analysis
	if not os.path.exists(path_analysis_pdb):
		os.makedirs(path_analysis_pdb)
	else:
		if len(os.listdir(path_analysis_pdb)) > 0:
			raise EnvironmentError("Analysis directory for Complex PDB contains files ")
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

def save_model_receptor(list_receptor_model_file):
	receptor_file = list_receptor_model_file[0]
	model_file = list_receptor_model_file[1]
	full_path_for_save_complex = list_receptor_model_file[2]

	#Open file for writting the complex
	f_compl = open(full_path_for_save_complex, "w")
	#Insert lines of receptor
	for item in  receptor_file:
		f_compl.write(item)
	#Insert lines of model
	for item in model_file:
		f_compl.write(item)
	f_compl.close()



