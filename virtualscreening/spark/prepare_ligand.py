from pyspark import SparkContext, SparkConf
import ConfigParser as configparser
from subprocess import Popen, PIPE
from datetime import datetime
import ntpath
import os
import vina_utils

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

def prepare_for_creating_database(database_path_file, pdbqt_path_file):
	if len(vina_utils.get_files_pdbqt(pdbqt_path_file)) == 0:
		raise EnvironmentError("pdbqt files not found when tried to create ligand Database")

	path = os.path.split(database_path_file)[0]
	if not os.path.exists(path):
		os.makedirs(path)	
"""
	Builds the ligand database
"""
def build_compound_database(f_pdbqt):
	ligand_name = get_name_ligand(f_pdbqt)
	torsion_angle = int(get_number_torsion_angle(f_pdbqt))
	atom_number = int(get_number_atom(f_pdbqt))
	line = (ligand_name + str("\t")+ str(torsion_angle) + str("\t") + str(atom_number) + str("\n") )	
	return line
	
"""
	Builds the ligand database
	database_path_file means the name of database. Here, contains path and name of file
"""
def save_database(database_path_file, list_compound):
	line = ""
	f_database = open(database_path_file, 'w')
	line = str(";Ligand\t") + str("Torsion\t") + str("Atoms\t") + str("\n")
	f_database.write(line)
	for line in list_compound:
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


""" This function gets the name of pdbqt files
"""
def get_name_ligand_pdbqt(reference):	
	name =  os.path.basename(reference)
	if str(name).find(".pdb") > 0:
		name = str(name).replace(".pdb", ".pdbqt")	
	if str(name).find(".mol2") > 0:
		name = str(name).replace(".mol2", ".pdbqt")	
	return name	


""" This function converts mol2 or pdb files 
    to pdbqt files for ligand
"""
def prepare_ligand(lig_vina):
		fpdbqt_filename = get_name_ligand_pdbqt(lig_vina[3])
		fpdbqt = os.path.join(lig_vina[0],fpdbqt_filename)				
		process = Popen([lig_vina[1], lig_vina[2], '-l', lig_vina[3], '-v', '-o', fpdbqt, '-U', 'nphs_lps', '-A', 'hydrogens'], stdout=PIPE, stderr=PIPE)		
		stdout, stderr = process.communicate()	 	

def save_log(finish_time, start_time):
	log_file_name = 'prepare_ligand.log'
	current_path = os.getcwd()
	path_file = os.path.join(current_path, log_file_name)
	log_file = open(path_file, 'w')

	diff_time = finish_time - start_time
	msg = 'Starting ' + str(start_time) +'\n'
	log_file.write(msg)
	msg = 'Finishing ' + str(finish_time) +'\n'
	log_file.write(msg)
	msg = 'Time Execution (seconds): ' + str(diff_time.total_seconds()) +'\n'
	log_file.write(msg)


def main():
	
	sc = SparkContext()
	config = configparser.ConfigParser()
	config.read('config.ini')

	#Broadcast - global
	path_pdbqt     = config.get('DEFAULT', 'pdbqt_ligand_path')
	pythonsh       = config.get('VINA', 'pythonsh')
	script_ligand4 = config.get('VINA', 'script_ligand4')
	database_comp  = config.get('DEFAULT', 'ligand_database_path_file')

	#creating pdbqt path
	if not os.path.isdir(path_pdbqt):
		os.mkdir(path_pdbqt)

	start_time = datetime.now()

	#preparing compound list
	list_obj_lig_vina = []
	mol2_files = vina_utils.get_files_mol2(config.get('DEFAULT', 'mol2_path'))
	for fmol2 in mol2_files:
		obj_lig_vina = (path_pdbqt, pythonsh,script_ligand4, fmol2)
		list_obj_lig_vina.append(obj_lig_vina)

	molRDD = sc.parallelize(list_obj_lig_vina)	
	molRDD.foreach(prepare_ligand)

	# *** Preparation of compound list finished. Now, it is able to create the database

	#preparing enviroment for creating database
	prepare_for_creating_database(database_comp, path_pdbqt)

	#preparing pdbqt list
	list_obj_pdbqt = []
	pdbqt_files = vina_utils.get_files_pdbqt(path_pdbqt)
	for fpdbqt in pdbqt_files:
		list_obj_pdbqt.append(fpdbqt)

	#appling map and collect
	pdbqtRDD = sc.parallelize(list_obj_pdbqt)
	all_lines = pdbqtRDD.map(build_compound_database).collect()
	

	#creating database file
	save_database(database_comp, all_lines)

	finish_time = datetime.now()

	save_log(finish_time, start_time)

main()




