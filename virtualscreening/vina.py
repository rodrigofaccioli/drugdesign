""" 
    Routines for interacting with AutoDock Vina program
    These routines were developed by:
    Rodrigo Antonio Faccioli - rodrigo.faccioli@usp.br / rodrigo.faccioli@gmail.com  
    Leandro Oliveira Bortot  - leandro.bortot@usp.br / leandro.obt@gmail.com 
"""



import sys
import os
import ntpath
from subprocess import Popen, PIPE
import shutil
import mol2

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

""" In this function is build the pdbqt file name for ligand 
"""
def get_name_ligand_pdbqt(reference):	
	name =  mol2.get_molecule_name(reference)
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
		fpdbqt_filename = get_name_ligand_pdbqt(fmol2)
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

""" This function obtains the name of receptor
	based on file name
"""	
def get_name_receptor(receptor):
	path, filename = ntpath.split(receptor)
	name =  str(filename.split(".")[0])
	return name

""" This function obtains the name of ligand
	based on file name
"""	
def get_name_ligand(ligand):
	path, filename = ntpath.split(ligand)
	name =  str(filename.split(".")[0])
	return name


""" This function creates out file name
"""	
def get_name_out(receptor, ligand):
	return receptor+'_-_'+ligand+'.pdbqt'

""" This function creates log file name
"""	
def get_name_log(receptor, ligand):
	return receptor+'_-_'+ligand+'.log'

""" This function is executed the docking from 
    one receptor against all ligands
"""	
def run_docking(vina_program, vina_conf, receptor, path_ligand_pdbqt, path_struct, path_log):	
	name_receptor = get_name_receptor(receptor)
	
	all_ligands = get_files_pdbqt(path_ligand_pdbqt)
	for ligand in all_ligands:
		name_ligand = get_name_ligand(ligand)		
		f_out = os.path.join(path_struct,get_name_out(name_receptor, name_ligand))
		f_log = os.path.join(path_log, get_name_log(name_receptor, name_ligand)) 			
		process = Popen([vina_program, '--config', vina_conf, '--receptor', receptor, '--ligand', ligand, '--out', f_out, '--log', f_log], stdout=PIPE, stderr=PIPE)		
		stdout, stderr = process.communicate()