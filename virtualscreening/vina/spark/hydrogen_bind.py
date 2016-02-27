from pyspark import SparkContext, SparkConf
from pyspark.sql import SQLContext, Row	
import ConfigParser as configparser
import os
import sys
from subprocess import Popen, PIPE
from datetime import datetime
from pdbqt_io import get_atom_section_from_atom_list, save_text_file_from_list
from vina_utils import get_name_model_pdb, get_name_model_pdbqt, get_directory_pdbqt_analysis, get_files_pdbqt, get_name_receptor_pdbqt, get_separator_filename_mode, get_ligand_from_receptor_ligand_model, get_model_from_receptor_ligand_model
import ntpath

def get_line_number(input_file):
	with open(input_file) as foo:
		lines_num = len( foo.readlines() )
	return lines_num

def get_saving_files_with_lines(mypath):
	only_saving_file = []
	for root, dirs, files in os.walk(mypath):
		for file in files:
			if file.endswith(".saving"):
				f_path = os.path.join(root,file)
				if get_line_number(f_path) > 0:
					only_saving_file.append(f_path)			
	return only_saving_file	

""" This function puts all lines of saving file
into list
"""
def loading_from_files(my_file_saving):
	list_return = []
	f_file = open(my_file_saving,"r")
	for line in f_file:		
		list_return.append( line ) 
	os.remove(my_file_saving)	
	return list_return

def loading_from_all_lists(sc, all_saving_filesRDD, sqlCtx):	
	#Splited all_saving_filesRDD into list that each element is a column.
	all_saving_filesRDD = sc.parallelize(all_saving_filesRDD)
	all_saving_filesRDD = all_saving_filesRDD.map(lambda s : str(s).split()).map(lambda p : Row(lig=str(p[0]), acpDon=str(p[1]), res=str(p[2]).strip(), atm=str(p[3]), hbondValue=float(p[4]), receptor=get_name_receptor_pdbqt(str(p[5])), ligand=get_ligand_from_receptor_ligand_model(str(p[6])), model=get_model_from_receptor_ligand_model( get_name_model_pdbqt(str(p[6]))) ))
	
			
	hbond_table = sqlCtx.createDataFrame(all_saving_filesRDD)
	hbond_table.registerTempTable("hbond")

	returnRDD = sqlCtx.sql("SELECT * FROM hbond ")
	return returnRDD

def get_all_saving_file(mypath):
	only_saving_file = []
	for root, dirs, files in os.walk(mypath):
		for file in files:
			if file.endswith(".saving"):
				f_path = os.path.join(root,file)
				only_saving_file.append(f_path)
	return only_saving_file

def remove_all_saving_files(mypath):
	all_saving_file = get_all_saving_file(mypath)
	for f in all_saving_file:
		os.remove(f)

""" This function obtains atom number, 3D position 
and atom ID (last column of pdbqt file) from atom list.
"""
def get_lig_values_from_atom_list_2_hydrogen_bind(str_atom_list_from_pdbqt):
	list_return = []
	for atom_line in str_atom_list_from_pdbqt:
		splited_line = atom_line.split()
		atom_num = splited_line[1]
		if len(splited_line) == 13:
			chain = 1
		else:
			chain = 0
		p_x = splited_line[5+chain]
		p_y = splited_line[6+chain]
		p_z = splited_line[7+chain]
		atomID = str(splited_line[-1])#Get the last column
		item = (atom_num, p_x, p_y, p_z, atomID)
		list_return.append(item)
	return list_return

""" This function obtains atom number, 3D position 
and atom ID (last column of pdbqt file) from atom list.
"""
def get_receptor_values_from_atom_list_2_hydrogen_bind(str_atom_list_from_pdbqt):
	list_return = []
	for atom_line in str_atom_list_from_pdbqt:
		splited_line = atom_line.split()
		if len(splited_line) == 13:
			chain = 1
		else:
			chain = 0
		res_name = splited_line[3]
		res_num = splited_line[4+chain]
		atom_name = splited_line[2]

		p_x = splited_line[5+chain]
		p_y = splited_line[6+chain]
		p_z = splited_line[7+chain]
		atomID = str(splited_line[-1])#Get the last column
		item = (res_name, res_num, atom_name, p_x, p_y, p_z, atomID)
		list_return.append(item)
	return list_return

def save_all_bonds_file(path_analysis, cutoff, all_saving_filesRDD):	
	f_file = "hbonds_all_"+str(cutoff)
	f_file = os.path.join(path_analysis, f_file)
	f_hbond = open(f_file,"w")	
	all_saving_filesRDD_2_txt = all_saving_filesRDD.map(lambda p: p.receptor + "\t"+ p.ligand +"\t"+ str(p.model)+"\t" + p.lig +"\t"+ p.acpDon +"\t"+ p.res +"\t"+ p.atm +"\t"+ str("{:.2f}".format(p.hbondValue)) +"\n")	
	for item in all_saving_filesRDD_2_txt.collect():
		f_hbond.write(item)
	f_hbond.close()

def save_all_no_bonds_file(path_analysis, path_saving_files, cutoff, sc):
	f_file = "NOT_hbonds_"+str(cutoff)
	f_file = os.path.join(path_analysis, f_file)
	f_no_hbond = open(f_file,"w")
	all_saving_file = get_all_saving_file(path_saving_files)

	RDD = sc.parallelize(all_saving_file).flatMap(lambda f : get_name_model_pdbqt(f)+"\n").collect()

	for line in  RDD:		
		f_no_hbond.write(line)
	f_no_hbond.close()

def get_hbonds_number_pose_constraints(constraints_file, path_analysis, sc, all_saving_filesRDD, sqlCtx):
	if os.path.isfile(constraints_file):		 
		f_constraints = sc.textFile(constraints_file)
		residues_consRDD = f_constraints.map(lambda l : str(str(l).split(";")[0]).strip() ).map(lambda p : Row(res=p))
		
		residues_cons_table = sqlCtx.createDataFrame(residues_consRDD)
		residues_cons_table.registerTempTable("residues_cons")

		number_pose_cons = sqlCtx.sql('SELECT hbond.receptor, hbond.ligand, hbond.model, count(hbond.res) as numPose FROM hbond INNER JOIN residues_cons ON residues_cons.res = hbond.res GROUP BY hbond.receptor, hbond.ligand, hbond.model')
		return number_pose_cons
	else:
		mensage = "This file was NOT found: \n"+constraints_file
		raise Exception(mensage)

def save_number_pose_constraints(path_analysis, cutoff, number_pose_consRDD):
	f_file = "hbonds_number_pose_constraints_"+str(cutoff)
	f_file = os.path.join(path_analysis, f_file)
	f_poses = open(f_file,"w")	
	all_poses_2_txt = number_pose_consRDD.map(lambda p: p.receptor +"\t"+ p.ligand +"\t"+ str(p.model) +"\t" + str(p.numPose) +"\n")	
	for item in all_poses_2_txt.collect():
		f_poses.write(item)
	f_poses.close()

def get_hbonds_number_pose(sqlCtx):
	
	number_pose = sqlCtx.sql('SELECT hbond.receptor, hbond.ligand, hbond.model, count(hbond.res) as numPose FROM hbond GROUP BY hbond.receptor, hbond.ligand, hbond.model')	
	return number_pose

def save_number_pose(path_analysis, cutoff, number_poseRDD):
	f_file = "hbonds_number_pose_"+str(cutoff)
	f_file = os.path.join(path_analysis, f_file)
	f_poses = open(f_file,"w")	
	all_poses_2_txt = number_poseRDD.map(lambda p: p.receptor +"\t"+ p.ligand +"\t"+ str(p.model) +"\t"+str(p.numPose) +"\n")	
	for item in all_poses_2_txt.collect():
		f_poses.write(item)
	f_poses.close()

def get_hbonds_number_ligand(sc, number_poseRDD, sqlCtx):
	only_ligandRDD = number_poseRDD.map(lambda p: Row(receptor=str(p[0]), ligand=str(p[1]), model=int(p[2]), numLigH=int(p[3]) ) )
	
	df = sqlCtx.createDataFrame(only_ligandRDD)
	df.registerTempTable("onlyligand")	
	
	number_ligand = sqlCtx.sql('SELECT onlyligand.receptor, onlyligand.ligand, SUM(numLigH) as num FROM onlyligand GROUP BY onlyligand.receptor, onlyligand.ligand')
	return number_ligand	

def save_number_ligand(path_analysis, cutoff, number_ligandRDD):
	f_file = "hbonds_number_ligand_"+str(cutoff)
	f_file = os.path.join(path_analysis, f_file)
	f_ligand = open(f_file,"w")	
	all_ligand_2_txt = number_ligandRDD.map(lambda p: p.receptor +"\t"+ p.ligand +"\t"+ str(p.num) +"\n")
	for item in all_ligand_2_txt.collect():
		f_ligand.write(item)
	f_ligand.close()


def save_vs_hydrogen_bind_log(finish_time, start_time):
	log_file_name = 'vs_hydrogen_bind.log'
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

def save_all_bonds_file_with_mensage(path_analysis, cutoff):
	f_file = "hbonds_all_"+str(cutoff)
	f_file = os.path.join(path_analysis, f_file)
	f_hbond = open(f_file,"w")	
	f_hbond.write("NO hydrogen bind")
	f_hbond.close()

def main():

	sc = SparkContext()
	sqlCtx = SQLContext(sc)

	config = configparser.ConfigParser()
	config.read('config.ini')

	#Path for drugdesign project
	path_spark_drugdesign = config.get('DRUGDESIGN', 'path_spark_drugdesign')
	#Detect interactions program
	detect_interactions_program = config.get('DRUGDESIGN', 'detect_interactions_program') 	
	#Path where all pdb receptor are
	path_receptor_pdbqt = config.get('DEFAULT', 'pdbqt_receptor_path')
	#Path that contains all files for analysis
	path_analysis = config.get('DEFAULT', 'path_analysis') 
	#Path of pdbqt model
	path_analysis_pdbqt_model = get_directory_pdbqt_analysis(path_analysis)

	#Getting parameters
	# cutoff for hydrogen bind
	cutoff = float(sys.argv[1])

	# checking for constraints file
	constraints_file = ""
	if len(sys.argv) >= 3:
		if sys.argv[2] != "":
			constraints_file = str(sys.argv[2])

	#Adding Python Source file
	sc.addPyFile(os.path.join(path_spark_drugdesign,"vina_utils.py"))
	sc.addPyFile(os.path.join(path_spark_drugdesign,"pdbqt_io.py"))

	start_time = datetime.now()

	#broadcast
	path_analysis_pdbqt_model_b = sc.broadcast(path_analysis_pdbqt_model)
	detect_interactions_program_b = sc.broadcast(detect_interactions_program)
	cutoff_b = sc.broadcast(cutoff)
#******************* start function ************************************************
	def get_hydrogen_bind(ligand_pdbqt):

		#getting base name
		base_name = get_name_model_pdb(ligand_pdbqt)
		#temporary_lig_no
		temporary_lig_no = base_name+"_temporary_lig_no"
		list_param = ["O", "N"]
		list_atom_pdbqt = get_atom_section_from_atom_list(ligand_pdbqt, list_param)	
		list_ref = get_lig_values_from_atom_list_2_hydrogen_bind(list_atom_pdbqt)
		path_file_lig_no = os.path.join(path_analysis_pdbqt_model_b.value, temporary_lig_no)
		save_text_file_from_list(path_file_lig_no, list_ref)
		total_lig_no = int(get_line_number(path_file_lig_no)) 

		#temporary_lig_h
		temporary_lig_h = base_name+"_temporary_lig_h"
		list_param = ["HD", "HS"]
		list_atom_pdbqt = get_atom_section_from_atom_list(ligand_pdbqt, list_param)	
		list_ref = get_lig_values_from_atom_list_2_hydrogen_bind(list_atom_pdbqt)
		path_file_lig_h = os.path.join(path_analysis_pdbqt_model_b.value, temporary_lig_h)
		save_text_file_from_list(path_file_lig_h, list_ref)
		total_lig_h = int(get_line_number(path_file_lig_h)) 

		#temporary_rec_no
		temporary_rec_no = base_name+"_temporary_rec_no"
		list_param = ["O", "N"]
		list_atom_pdbqt = get_atom_section_from_atom_list(receptor_b.value, list_param)	
		list_ref = get_receptor_values_from_atom_list_2_hydrogen_bind(list_atom_pdbqt)
		path_file_rec_no = os.path.join(path_analysis_pdbqt_model_b.value, temporary_rec_no)
		save_text_file_from_list(path_file_rec_no, list_ref)
		total_rec_no = int(get_line_number(path_file_rec_no)) 

		#temporary_rec_h
		temporary_rec_h = base_name+"_temporary_rec_h"
		list_param = ["HD", "HS"]
		list_atom_pdbqt = get_atom_section_from_atom_list(receptor_b.value, list_param)	
		list_ref = get_receptor_values_from_atom_list_2_hydrogen_bind(list_atom_pdbqt)
		path_file_rec_h = os.path.join(path_analysis_pdbqt_model_b.value, temporary_rec_h)
		save_text_file_from_list(path_file_rec_h, list_ref)
		total_rec_h = int(get_line_number(path_file_rec_h)) 
		
		#preparing file for saving	
		file_for_saving = base_name+".saving"
		path_file_for_saving = os.path.join(path_analysis_pdbqt_model_b.value, file_for_saving)		
		if total_lig_h > 0 and total_lig_no > 0:		
			process = Popen( [detect_interactions_program_b.value, receptor_b.value, str(total_rec_no), str(total_rec_h), ligand_pdbqt, str(total_lig_no), str(total_lig_h), str(cutoff_b.value),path_file_for_saving,path_file_rec_no, path_file_lig_no, path_file_rec_h, path_file_lig_h ], stdout=PIPE, stderr=PIPE)
			stdout, stderr = process.communicate()

		#os.remove(path_file_rec_no)
		#os.remove(path_file_lig_no)
		#os.remove(path_file_rec_h)
		#os.remove(path_file_lig_h)

#******************* finish function ************************************************

	#Getting all receptores
	all_receptores = get_files_pdbqt(path_receptor_pdbqt)

	#Getting all pdbqt models
	all_pdbqt_models = get_files_pdbqt(path_analysis_pdbqt_model)
	all_pdbqt_modelsRDD = sc.parallelize(all_pdbqt_models)

	for receptor in all_receptores:
		receptor_b = sc.broadcast(receptor)
		base_name_receptor = get_name_receptor_pdbqt(receptor)
		models_by_receptorRDD = all_pdbqt_modelsRDD.filter(lambda m : base_name_receptor in m).collect()
		models_by_receptorRDD = sc.parallelize(models_by_receptorRDD)
		models_by_receptorRDD.foreach(get_hydrogen_bind)

	#Getting all saving files that have lines > 0
	all_saving_files = get_saving_files_with_lines(path_analysis_pdbqt_model)

	if len(all_saving_files) > 0:

		all_saving_filesRDD = sc.parallelize(all_saving_files)

		#loading from files
		all_saving_filesRDD = all_saving_filesRDD.flatMap(loading_from_files).collect()

		#loading all values from list
		all_saving_filesRDD = loading_from_all_lists(sc, all_saving_filesRDD, sqlCtx)
		all_saving_filesRDD.cache()

		#saving all_bonds_file	
		save_all_bonds_file(path_analysis, cutoff, all_saving_filesRDD)

		#saving models of ligands which no have hydrogen bind
		save_all_no_bonds_file(path_analysis, path_analysis_pdbqt_model, cutoff, sc)

		#removing remaing saving files (They no lines have)
		remove_all_saving_files(path_analysis_pdbqt_model)	

		#number hydrogen binds of poses by constraints
		if constraints_file != "":
			number_pose_consRDD = get_hbonds_number_pose_constraints(constraints_file, path_analysis, sc, all_saving_filesRDD, sqlCtx)
			save_number_pose_constraints(path_analysis, cutoff, number_pose_consRDD)

		#number hydrogen binds of poses
		number_poseRDD = get_hbonds_number_pose(sqlCtx)
		number_poseRDD.cache()
		save_number_pose(path_analysis, cutoff, number_poseRDD)

		#number hydrogen binds of ligands
		number_ligandRDD = get_hbonds_number_ligand(sc, number_poseRDD, sqlCtx)
		save_number_ligand(path_analysis, cutoff, number_ligandRDD)

	else:
		save_all_bonds_file_with_mensage(path_analysis, cutoff)
		#removing remaing saving files (They no lines have)
		remove_all_saving_files(path_analysis_pdbqt_model)	


	finish_time = datetime.now()

	save_vs_hydrogen_bind_log(finish_time, start_time)

main()