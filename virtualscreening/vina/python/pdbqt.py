""" 
    Routines for manipulating pdbqt files
    These routines were developed by:
    Rodrigo Antonio Faccioli - rodrigo.faccioli@usp.br / rodrigo.faccioli@gmail.com  
    Leandro Oliveira Bortot  - leandro.bortot@usp.br / leandro.obt@gmail.com 
"""

import os
import ntpath
import operator

""" This function obtains the name of ligand
	based on file name
"""	
def get_name_ligand(ligand):
	path, filename = ntpath.split(ligand)
	name =  str(filename.split(".")[0])
	return name

"""
	Returns the number of torsion angles 
	
	pathfilename contains path and file name in which
	wants to know the number of torsion angles.
	In pdbqt files, the number of torsion angles is 
	REMARK  2 active torsions	
"""
def get_number_torsion_angle(pathfilename):
	number = int(0)
	fpdbqt = open(pathfilename, "r")
	for line in fpdbqt:
		if (line.find("active torsions") > 0):
			s = str(line).split()
			number = int(s[1])
			break
	fpdbqt.close()
	return number

def sort_dictionary(dic_docking):
	return sorted(dic_docking.items(), key=operator.itemgetter(1), reverse=True)

def get_files_pdbqt_no_sort(mypath):
	only_pdbqt_file = []
	for root, dirs, files in os.walk(mypath):
		for file in files:
			if file.endswith(".pdbqt"):
				f_path = os.path.join(root,file)
				only_pdbqt_file.append(f_path)
	return only_pdbqt_file

"""
	Returns the number of atoms 
	
	pathfilename contains path and file name in which
	wants to know the number of atoms.
	In pdbqt files, the number of atom is obatined reading
	all lines of files looking by ATOM. The number of atom 
	is the greater number.
"""
def get_number_atom(pathfilename):	
	number = int(0)
	fpdbqt = open(pathfilename, "r")
	for line in fpdbqt:
		if (line.find("ATOM") > -1):			
			s = str(line).split()
			number = int(s[1])			
	fpdbqt.close()
	return number