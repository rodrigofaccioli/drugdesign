#! /usr/bin/env python
""" 
    Routines to extract compounds from ZINC Database: http://zinc.docking.org/
    These routines were developed by:
    Rodrigo Antonio Faccioli - rodrigo.faccioli@usp.br / rodrigo.faccioli@gmail.com  
    Leandro Oliveira Bortot  - leandro.bortot@usp.br / leandro.obt@gmail.com 
"""

import ConfigParser as configparser
import os
import shutil 
import gzip

def number_files_of_molecule(molecule_name, path_save_mol2):
	"""
    Return the number of files at path_save_mol2 that contain molecule_name in filename
    Example:
        >>> number = number_files_of_molecule(molecule_name, path_save_mol2)
    @param molecule_name: main name of molecule
    @type molecule_name: string
    @param path_save_mol2: path of mol2 files will be saved
    @type path_save_mol2: string    
    @return: the number of files that have molecule_name in their names
    @rtype: int    
	"""
	number = 0
	for root, dirs, files in os.walk(path_save_mol2):
		for f in files:			
			if f.endswith(".mol2"): 				
				if str(f).find(molecule_name) >=0:
					number = number + 1		
	return number

def finish_current_molecule(molecule_name, path_save_mol2, temp_file_name_full):
	"""
    Last procedures for current molecule  
    Example:
        >>> finish_current_molecule(molecule_name, path_save_mol2, temp_file_name_full)
    @param molecule_name: main name of molecule
    @type molecule_name: string
    @param path_save_mol2: path of mol2 files will be saved
    @type path_save_mol2: string    
    @param temp_file_name_full: full path of temp file
    @type temp_file_name_full: string
	"""		
	#preparing name of mol2 file
	# Checking filenames of molecules based on molecule_name
	# Because of isomers, it is necessary to check how many files of 
	# molecule_name there is in path_save_mol2		
	mol_name_aux = ''
	number_files = number_files_of_molecule(molecule_name, path_save_mol2)	
	if  number_files > 0:
		if number_files == 1:
			#means that there is only one molecule.
			#So it must be renamed with prefix _1
			#number_files will be assigned to 2, because
			# the current molecule will be second molecule			
			before_molecule = molecule_name+'.mol2'
			before_molecule_mol2 = os.path.join(path_save_mol2, before_molecule)
			new_molecule = molecule_name+'_1'+'.mol2'
			new_molecule_mol2 = os.path.join(path_save_mol2, new_molecule)			
			shutil.move(before_molecule_mol2, new_molecule_mol2)
		number_files = number_files + 1			
		mol_name_aux = molecule_name+'_'+str(number_files)
	else:
		mol_name_aux = molecule_name
	mol2_file_name = mol_name_aux+'.mol2'
	mol2_file_name_full = os.path.join(path_save_mol2, mol2_file_name)
	#creating mol2 file - moving temp file to mol2_file_name_full			
	shutil.move(temp_file_name_full, mol2_file_name_full)

def split_molecules_from_mol2(pathfilename, path_save_mol2):
	"""
    Split molecules from mol2 file
    Example:
        >>> split_molecules_from_mol2(pathfilename, path_save_mol2)
    @param pathfilename: full path of file that contains all molecules
    @type pathfilename: string
    @param path_save_mol2: path of mol2 files will be saved
    @type path_save_mol2: string    
	"""	
	line_name = True
	new_molecule = False	
	temp_file_name = 'temp.temp'
	line_aux = ""

	#open full mol2 file
	fmol2_all = open(pathfilename, "r")
	##open temp file for first molecule
	temp_file_name_full = os.path.join(path_save_mol2, temp_file_name)
	fmol2_temp = open(temp_file_name_full, "w")
	#Obtain first line from full mol2 file
	line = fmol2_all.readline()
	fmol2_temp.write(line)
	for line in fmol2_all:
		if line.find("@<TRIPOS>MOLECULE") < 0:
			#get the molecule name
			if line_name == True:
				molecule_name = str(line).strip()
				line_name = False			
			fmol2_temp.write(line)
		else: # found @<TRIPOS>MOLECULE
			#close temp file
			fmol2_temp.close()
			#finishing the current molecule
			finish_current_molecule(molecule_name, path_save_mol2, temp_file_name_full)			
			#open temp file for new molecue
			temp_file_name_full = os.path.join(path_save_mol2, temp_file_name)
			fmol2_temp = open(temp_file_name_full, "w")
			#assign True to line_name
			line_name = True
			#assign line to temp file
			fmol2_temp.write(line)
	#close temp file
	fmol2_temp.close()
	#finishing the last molecule
	finish_current_molecule(molecule_name, path_save_mol2, temp_file_name_full)			

def get_files_gz(mypath):
	only_gz_file = []
	for root, dirs, files in os.walk(mypath):
		for file in files:
			if file.endswith(".gz"):
				f_path = os.path.join(root,file)
				only_gz_file.append(f_path)			
	return only_gz_file

def decompress_gz_files(gz_path):
	f_name = ""
	path_filename = ""
	gz_files = get_files_gz(gz_path)
	for f in gz_files:
		f_name = os.path.basename(f)
		f_name = str(f_name).replace(".gz",'')
		path_filename = os.path.join(gz_path, f_name)
		inF = gzip.open(f, 'rb')
		s = inF.read()
		inF.close()
		mol2_file = open(path_filename, 'w')
		mol2_file.write(s)
		mol2_file.close()

def get_files_mol2(mypath):
	only_mol2_file = []
	for root, dirs, files in os.walk(mypath):
		for file in files:
			if file.endswith(".mol2"):
				f_path = os.path.join(root,file)
				only_mol2_file.append(f_path)			
	return only_mol2_file		

def main():
	config = configparser.ConfigParser()
	config.read('config.ini')

	decompress_gz_files(config.get('ZINCDB', 'path_downloaded'))
	mol2_files = get_files_mol2(config.get('ZINCDB', 'path_downloaded'))
	for f_mol2 in mol2_files:
		split_molecules_from_mol2(f_mol2, config.get('DEFAULT', 'mol2_path') )

main()