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
