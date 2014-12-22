""" 
    Routines to preparation of enviroment for virtual screening execution
    These routines were developed by:
    Rodrigo Antonio Faccioli - rodrigo.faccioli@usp.br / rodrigo.faccioli@gmail.com  
    Leandro Oliveira Bortot  - leandro.bortot@usp.br / leandro.obt@gmail.com 
"""

import os

def get_files_pdbqt(mypath):
	only_pdbqt_file = []
	for root, dirs, files in os.walk(mypath):
		for file in files:
			if file.endswith(".pdbqt"):
				f_path = os.path.join(root,file)
				only_pdbqt_file.append(f_path)			
	return only_pdbqt_file

def preparing_vs(path_save_structure, path_save_log, path_analysis, pdbqt_ligand_path, pdbqt_receptor_path):

	#Checking path_save_structure
	if not os.path.exists(path_save_structure):
		os.makedirs(path_save_structure)
	else:
		if len(os.listdir(path_save_structure)) > 0:
			raise EnvironmentError("Structure directory contains files ")

	#Checking path_save_log
	if not os.path.exists(path_save_log):
		os.makedirs(path_save_log)
	else:
		if len(os.listdir(path_save_log)) > 0:
			raise EnvironmentError("Log directory contains files ")

	#Checking path_analysis
	if not os.path.exists(path_analysis):
		os.makedirs(path_analysis)
	else:
		if len(os.listdir(path_analysis)) > 0:
			raise EnvironmentError("Analysis directory contains files ")

	#Checking pdbqt files of ligand
	if len(get_files_pdbqt(pdbqt_ligand_path)) == 0:
		raise EnvironmentError("pdbqt files of ligands not found. Remember, prepare_ligand.py must be executed ")

	#Checking pdbqt files of receptor
	if len(get_files_pdbqt(pdbqt_receptor_path)) == 0:
		raise EnvironmentError("pdbqt files of receptors not found. Remember, prepare_receptor.py must be executed ")