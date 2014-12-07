# functions to call Auto Vina Docking program


import sys
import os
import ntpath
from subprocess import Popen, PIPE
import shutil

""" This function obtains all pdbqt files 
    in mypath 
"""
def get_files_pdbqt(mypath):
	only_pdbqt_file = []
	for root, dirs, files in os.walk(mypath):
		for file in files:
			if file.endswith(".pdbqt"):
				f_path = os.path.join(root,file)
				only_pdbqt_file.append(f_path)			
	return only_pdbqt_file


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

""" This function obtains all mol2 files 
    in mypath 
"""
def get_files_mol2(mypath):
	only_mol2_file = []
	for root, dirs, files in os.walk(mypath):
		for file in files:
			if file.endswith(".mol2"):
				f_path = os.path.join(root,file)
				only_mol2_file.append(f_path)			
	return only_mol2_file

""" This function is created the pdbqt file get_name_pdbqt
     based on fmol2 file name 
"""
def get_name_pdbqt(reference):
	path, filename = ntpath.split(reference)
	name =  str(filename.split(".")[0]) #remove either .mol2 or .pdb
	fpdbqt = name+".pdbqt"
	return fpdbqt


""" This function converts mol2 files 
    to pdbqt files 
"""
def prepare_ligand(path_mol2, path_pdbqt, pythonsh, script_ligand4):

	if not os.path.isdir(path_pdbqt):
		os.mkdir(path_pdbqt)

	mol2_files = get_files_mol2(path_mol2)
	for fmol2 in mol2_files:
		fpdbqt_filename = get_name_pdbqt(fmol2)
		fpdbqt = os.path.join(path_pdbqt,fpdbqt_filename)
		process = Popen([pythonsh, script_ligand4, '-l', fmol2, '-v', '-o', fpdbqt, '-U', 'nphs_lps', '-A', 'hydrogens'], stdout=PIPE, stderr=PIPE)		
		stdout, stderr = process.communicate()	 	
		#command = pythonsh + " " + script_ligand4 + " "+ '-l' + " "+ fmol2 + " "+ '-v' + " "+ '-o'+ " "+ fpdbqt + "\n"
		#os.system(command)

""" This function converts pdb files 
    to pdbqt files 
"""	
def prepare_receptor(path_pdb, path_pdbqt, pythonsh, script_receptor4):
	if not os.path.isdir(path_pdbqt):
		os.mkdir(path_pdbqt)

	pdb_files = get_files_pdb(path_pdb)
	for fpdb in pdb_files:
		fpdbqt_filename = get_name_pdbqt(fpdb)
		fpdbqt = os.path.join(path_pdbqt,fpdbqt_filename)
		process = Popen([pythonsh, script_receptor4, '-r', fpdb, '-o', fpdbqt, '-v', '-U', 'nphs_lps', '-A', 'hydrogens'], stdout=PIPE, stderr=PIPE)		
		stdout, stderr = process.communicate()	 	

""" This function creates out file name
"""	
def get_name_out(index):
	return 'estrutura_'+str(index)+'.pdbqt'

""" This function creates log file name
"""	
def get_name_log(index):
	return 'log_'+str(index)+'.log'

""" This function is executed the docking from 
    one receptor against all ligands
"""	
def run_docking(vina_program, vina_conf, receptor, path_ligand_pdbqt, path_struct, path_log):
	index = int(-1)
	all_ligands = get_files_pdbqt(path_ligand_pdbqt)
	for ligand in all_ligands:
		index = index + 1
		f_out = os.path.join(path_struct,get_name_out(index))
		f_log = os.path.join(path_log, get_name_log(index)) 
		process = Popen([vina_program, '--config', vina_conf, '--receptor', receptor, '--ligand', ligand, '--out', f_out, '--log', f_log], stdout=PIPE, stderr=PIPE)		
		stdout, stderr = process.communicate()