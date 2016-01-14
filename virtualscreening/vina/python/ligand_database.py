#! /usr/bin/env python
""" 
    Routines for preparing the ligand database
    These routines were developed by:
    Rodrigo Antonio Faccioli - rodrigo.faccioli@usp.br / rodrigo.faccioli@gmail.com  
    Leandro Oliveira Bortot  - leandro.bortot@usp.br / leandro.obt@gmail.com 
"""

import os
import pdbqt

def prepare_for_creating_database(database_path_file, pdbqt_path_file):
	if len(pdbqt.get_files_pdbqt_no_sort(pdbqt_path_file)) == 0:
		raise EnvironmentError("pdbqt files not found when tried to create ligand Database")

	path = os.path.split(database_path_file)[0]
	if not os.path.exists(path):
		os.makedirs(path)	

"""
	Builds the ligand database
	
	database_path_file means the name of database. Here, contains path and name of file
	pdbqt_path_file means path where the pdbqt files are
"""
def build_database(database_path_file, pdbqt_path_file):
	line = ""
	f_database = open(database_path_file, 'w')
	line = str(";Ligand\t") + str("Torsion\t") + str("Atoms\t") + str("\n")
	f_database.write(line)
	all_pdbqt = pdbqt.get_files_pdbqt_no_sort(pdbqt_path_file)
	for f_pdbqt in all_pdbqt:
		ligand_name = pdbqt.get_name_ligand(f_pdbqt)
		torsion_angle = int(pdbqt.get_number_torsion_angle(f_pdbqt))
		atom_number = int(pdbqt.get_number_atom(f_pdbqt))
		line = ligand_name + str("\t")+ str(torsion_angle) + str("\t") + str(atom_number) + str("\n")
		f_database.write(line)

	f_database.close()

"""
	Returns torsion and atom number from ligand database
	
	compound_name means the name of compound (ligand) that
	wants to obtain the number of torsion angles and the number
	of atoms
"""		
def get_torsion_atom_from_database(compound_name, ligand_database):
	torsion = int(0)
	atom = int(0)
	database_file = open(ligand_database, 'r')
	for line in database_file:		
		if str(line).find(compound_name) > -1:
			splited_line = str(line).split()
			torsion = int(splited_line[1])
			atom = int(splited_line[2])
			break
	database_file.close()
	return torsion, atom
