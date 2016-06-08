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

""" This function obtains all pdb files 
    in mypath filtered by reference
"""
def get_files_pdb_filter(mypath, reference):
	only_pdb_file = []
	for root, dirs, files in os.walk(mypath):
		for file in files:
			if file.endswith(".pdb"):
				if file.find(reference) >  -1:
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
	return 'vs_sorted_energies_vina.dat'


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
	return path_analysis_pdbqt

""" This function obtains the name of 
path that saving pdbqt files for analysis 
"""
def get_directory_temp_analysis(path_analysis):
	path_analysis_temp = os.path.join(path_analysis,"temp")
	#Checking path_analysis
	if not os.path.exists(path_analysis_temp):
		os.makedirs(path_analysis_temp)
	return path_analysis_temp


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

def get_name_model_pdbqt(myfile):
	""" 
	This function obtains the name of myfile without filename extension
	"""	
	path, filename = ntpath.split(myfile)
	name =  str(filename.split(".")[0]) #remove .pdbqt
	return name

def get_name_receptor_pdb(myfile):
	""" 
	This function obtains the name of myfile without filename extension
	"""	
	path, filename = ntpath.split(myfile)
	name =  str(filename.split(".")[0]) #remove .pdb
	return name

def get_name_receptor_pdbqt(myfile):
	""" 
	This function obtains the name of myfile without filename extension
	"""	
	path, filename = ntpath.split(myfile)
	name =  str(filename.split(".")[0]) #remove .pdbqt
	return name

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


def get_ligand_from_receptor_ligand_model(receptor_ligand_model):
	""" 
	This function obtains the name of ligand based on receptor_ligand_model
	Example of input: compl_ns3pro_dm_0_-_NuBBE_485_obabel_3D+----+20
	"""	
	separator_model = get_separator_filename_mode()
	separator_receptor = "_-_"
	string_ref = receptor_ligand_model
	
	s = string_ref.split(separator_receptor) #Removing receptor
	s = str(s[1]).split(separator_model) #Removing model
	ligand_name = str(s[0]) #geting name of ligand
	return ligand_name

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
	model = int(s[1]) #geting model
	return model


""" This function obtains the name of 
path that saving pdbqt files for analysis 
"""
def get_directory_pdb_analysis(path_analysis):
	path_analysis_pdb = os.path.join(path_analysis,"pdb_model")
	#Checking path_analysis
	if not os.path.exists(path_analysis_pdb):
		os.makedirs(path_analysis_pdb)
	return path_analysis_pdb

""" This function obtains the name of 
path that saving pdbqt files for analysis 
"""
def get_directory_complex_pdb_analysis(path_analysis):
	path_analysis_pdb = os.path.join(path_analysis,"pdb_complex")
	#Checking path_analysis
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

