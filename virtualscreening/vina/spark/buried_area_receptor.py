from pyspark import SparkContext, SparkConf, SparkFiles
from pyspark.sql import SQLContext, Row
import ConfigParser as configparser
from subprocess import Popen, PIPE
from datetime import datetime
from vina_utils import get_directory_complex_pdb_analysis, get_files_pdb, get_name_model_pdb, get_ligand_from_receptor_ligand_model, get_separator_filename_mode, get_directory_pdb_analysis, loading_pdb_2_list, get_name_receptor_pdb, get_files_pdb_filter, get_receptor_from_receptor_ligand_model, get_model_from_receptor_ligand_model
import os, sys
from os_utils import preparing_path
from gromacs_utils import get_value_from_xvg_sasa
from pdb_io import replace_chain_atom_line
import shutil

def sorting_by_buried_lig_rec(sc, buried_areaRDD):
	sqlCtx = SQLContext(sc)
	buried_areaRDD = sc.parallelize(buried_areaRDD)
	buried_areaRDD = buried_areaRDD.map(lambda p: Row(pose=str(p[0]), buried_lig_rec=float(p[1]) ) )
	buried_area_table = sqlCtx.createDataFrame(buried_areaRDD)
	buried_area_table.registerTempTable("buried_area_receptor_summary")

	buried_area_sorted_by_lig_rec = sqlCtx.sql("SELECT pose, buried_lig_rec FROM buried_area_receptor_summary ORDER BY buried_lig_rec DESC")
	return buried_area_sorted_by_lig_rec

def sorting_buried_area_receptor(sc, buried_areaRDD):
	sqlCtx = SQLContext(sc)
	buried_areaRDD = sc.parallelize(buried_areaRDD)
	buried_areaRDD = buried_areaRDD.map(lambda p: Row(pose=str(p[0]), buried_area_rec=float(p[1]) ) )
	buried_area_table = sqlCtx.createDataFrame(buried_areaRDD)
	buried_area_table.registerTempTable("buried_area_receptor")

	buried_area_sorted_by_recep = sqlCtx.sql("SELECT * FROM buried_area_receptor ORDER BY buried_area_rec ")
	return buried_area_sorted_by_recep

def save_buried_area_receptor_sort(path_file_buried_area, buried_area):
	f_buried_area = open(path_file_buried_area,"w")
	line = "# buried_area_receptor[nm2]\tpose"+"\n"
	f_buried_area.write(line)
	for area in buried_area:
		pose = str(str(area[0]).replace("compl_", " ")).strip()
		buried_lig_rec = "{:.4f}".format(area[1])
		line = str(buried_lig_rec)+"\t"+str(pose)+"\n"
		f_buried_area.write(line)
	f_buried_area.close()


def sorting_buried_area_all_residues(sc, buried_areaRDD):
	sqlCtx = SQLContext(sc)
	buried_areaRDD = sc.parallelize(buried_areaRDD)
	#buried_areaRDD = buried_areaRDD.map(lambda p: Row(receptor=str(p[0]), ligand=str(p[1]), model=int(p[2]), res=str(p[3]), res_buried_area=float(p[4]), res_buried_area_perc=float(p[5]) ) )
	buried_areaRDD = buried_areaRDD.map(lambda p: Row(res=str(p[0]), res_buried_area=float(p[1]), res_buried_area_perc=float(p[2]), pose=str(p[3]) ) )
	buried_area_table = sqlCtx.createDataFrame(buried_areaRDD)
	buried_area_table.registerTempTable("buried_area_recep")

	buried_area_sorted_by_lig_rec_perc = sqlCtx.sql("SELECT * FROM buried_area_recep ORDER BY pose, res_buried_area ") #receptor,ligand,model,
	return buried_area_sorted_by_lig_rec_perc

def save_buried_area_recep(path_file_buried_area, buried_area_sorted_by_res_buried_area_perc):
	f_buried_area = open(path_file_buried_area,"w")
	line = "# residue\tburied_area_residue[nm2]\tresidue_sasa_buried[%]\tpose"+"\n"
	f_buried_area.write(line)
	for area in buried_area_sorted_by_res_buried_area_perc:
		#receptor = area[0]
		#ligand = area[1]
		#model = area[2]
		res = str(area[0]).replace("_","-")
		res_buried_area = "{:.4f}".format(area[1])
		res_buried_area_perc = "{:.4f}".format(area[2])
		pose = area[3]
		#line = receptor+"\t"+ligand+"\t"+str(model)+"\t"+str(res)+"\t"+str(res_buried_area)+"\t"+str(res_buried_area_perc)+"\n"
		line = str(res)+"\t"+str(res_buried_area)+"\t"+str(res_buried_area_perc)+"\t"+str(pose)+"\n"
		f_buried_area.write(line)
	f_buried_area.close()


def save_receptor_buried_area_receptor(path_file_buried_area, buried_area):
	f_buried_area = open(path_file_buried_area,"w")
	for area in buried_area:
		for t in area:
			#splited_line = t[0].split("_-_")
			#aux_recep = splited_line[0]
			#aux_lig = str(splited_line[1])
			#preparing receptor
			#receptor = str(str(aux_recep).replace("compl_", " ")).strip()
			#preparing ligand
			#splited_aux_lig = str(aux_lig).split(get_separator_filename_mode())
			#ligand = splited_aux_lig[0]
			#model = splited_aux_lig[1]
			pose = str(str(t[0]).replace("compl_", " ")).strip()
			res = t[1]
			res_buried_area = "{:.4f}".format(t[2])
			res_buried_area_perc = "{:.4f}".format(t[3])
			#line = receptor+"\t"+ligand+"\t"+model+"\t"+res+"\t"+str(res_buried_area)+"\t"+str(res_buried_area_perc)+"\n"
			line = str(res)+"\t"+str(res_buried_area)+"\t"+str(res_buried_area_perc)+"\t"+str(pose)+"\n"
			f_buried_area.write(line)
	f_buried_area.close()

def loading_lines_from_recepArea_files(line):
	line_splited = str(line).split()
	#line_ret = ( str(line_splited[0]), str(line_splited[1]), int(line_splited[2]), str(line_splited[3]), float(line_splited[4]), float(line_splited[5]) )
	line_ret = ( str(line_splited[0]), float(line_splited[1]), float(line_splited[2]), str(line_splited[3]) )
	return line_ret

def loading_lines_from_outAreaRecep_files(line):
	line_splited = str(line).split()
	line_ret = ( str(line_splited[0]), float(line_splited[1]) )
	return line_ret

def get_files_recepArea(mypath):
	only_mol2_file = []
	for root, dirs, files in os.walk(mypath):
		for file in files:
			if file.endswith(".recepArea"):
				f_path = os.path.join(root,file)
				only_mol2_file.append(f_path)
	return only_mol2_file


def get_files_outAreaRecep(mypath):
	only_mol2_file = []
	for root, dirs, files in os.walk(mypath):
		for file in files:
			if file.endswith(".outAreaRecep"):
				f_path = os.path.join(root,file)
				only_mol2_file.append(f_path)
	return only_mol2_file


def get_residues_receptor_from_ndx_files(f_name_ndx):
	list_res = []
	f_ndx = open(f_name_ndx,"r")
	for line in f_ndx:
		if line.find("rec_") > -1:
			line = line.translate(None, "]")
			line = line.translate(None, "[")
			line = line.strip()
			res = line.split("_")
			res = res[1]+"_"+res[2]
			list_res.append(res)
	f_ndx.close()
	return list_res

def save_log(finish_time, start_time):
	log_file_name = 'vs_buried_areas_receptor.log'
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

	config = configparser.ConfigParser()
	config.read('config.ini')

	#Path for Gromacs project
	gromacs_path = preparing_path(config.get('DRUGDESIGN', 'gromacs_path'))
	#Path where PDB ligand are - They are NOT participated in docking
	pdb_ligand_path = config.get('DEFAULT', 'pdb_ligand_path')
	#Path that contains all files for analysis
	path_analysis = config.get('DEFAULT', 'path_analysis')
	#Path where all pdb receptor are
	path_receptor_pdb = config.get('DEFAULT', 'pdb_path')
	#Path for saving pdb files of models generated by VS
	path_analysis_pdb = get_directory_pdb_analysis(path_analysis)

	# Create SPARK config
	maxResultSize = str(config.get('SPARK', 'maxResultSize'))
	conf = (SparkConf().set("spark.driver.maxResultSize", maxResultSize))

	# Create context
	sc = SparkContext(conf=conf)

	#Adding Python Source file
	#Path for drugdesign project
	path_spark_drugdesign = config.get('DRUGDESIGN', 'path_spark_drugdesign')
	sc.addPyFile(os.path.join(path_spark_drugdesign,"vina_utils.py"))
	sc.addPyFile(os.path.join(path_spark_drugdesign,"os_utils.py"))
	sc.addPyFile(os.path.join(path_spark_drugdesign,"gromacs_utils.py"))
	sc.addPyFile(os.path.join(path_spark_drugdesign,"pdb_io.py"))
	sc.addPyFile(os.path.join(path_spark_drugdesign,"json_utils.py"))

	#Adding bash scripts
	sc.addFile(os.path.join(path_spark_drugdesign,"make_ndx_buried_area_receptor.sh"))
	sc.addFile(os.path.join(path_spark_drugdesign,"make_ndx_buried_area_receptor_res.sh"))

	#Parameters form command line
	#Indicates probe. Example: 0.14
	#probe = float(sys.argv[1])
	#Indicates ndots. Example: 24
	#ndots = int(sys.argv[2])

	#Broadcast
	path_analysis_pdb_complex_b = sc.broadcast(path_analysis_pdb)
	gromacs_path = sc.broadcast(gromacs_path)
	pdb_ligand_path = sc.broadcast(pdb_ligand_path)
	#probe = sc.broadcast(probe)
	#ndots = sc.broadcast(ndots)

	start_time = datetime.now()

	os.environ["GMX_MAXBACKUP"]="-1"

	#Loading all PDB receptor files into memory
	list_all_pdb_receptor_files_path = []
	all_receptor_for_complex = get_files_pdb(path_receptor_pdb)
	for receptor in all_receptor_for_complex:
		list_all_pdb_receptor_files_path.append(loading_pdb_2_list(receptor))

	for pdb_receptor_files in list_all_pdb_receptor_files_path:
		#Getting receptor name by fully path
		base_file_name_receptor = get_name_receptor_pdb(str(pdb_receptor_files[0]))
		#PDB file loaded into memory is sent by broadcast
		pdb_file_receptor = pdb_receptor_files[1]
		pdb_file_receptor = sc.broadcast(pdb_file_receptor)
		#Loading PDB model files based on receptor into memory
		base_file_name_receptor_for_filter = base_file_name_receptor+"_-_"
		all_model_for_complex = get_files_pdb_filter(path_analysis_pdb,base_file_name_receptor_for_filter)
		all_model_for_complexRDD = sc.parallelize(all_model_for_complex)
		all_model_filesRDD = all_model_for_complexRDD.map(loading_pdb_2_list).collect()

# ********** Starting function **********************************************************
		def save_model_receptor(list_receptor_model_file):
			receptor_file = pdb_file_receptor.value #Obtained from broadcast
			model_file = list_receptor_model_file[0]
			full_path_for_save_complex = list_receptor_model_file[1]
			#Open file for writting the complex
			f_compl = open(full_path_for_save_complex, "w")
			#Insert lines of receptor
			for item in  receptor_file:
				f_compl.write(item)
			#Insert lines of model and insert Z chain
			for item in model_file:
				item = replace_chain_atom_line(item,"d","z")
				f_compl.write(item)
			f_compl.close()
# ********** Finish function **********************************************************

# ********** Starting function **********************************************************
		def compute_buried_area_all_residues_and_receptor_area(pdb_complex):
			chZ = "chZ"
			res_buried_area_perc = -1
			res_buried_area = -1
			buried_receptor_system = -1
			buried_receptor_res = -1
			base_name = get_name_model_pdb(pdb_complex)
			ligand_name = get_ligand_from_receptor_ligand_model(base_name)
			receptor_name = get_receptor_from_receptor_ligand_model(base_name)
			pose = get_model_from_receptor_ligand_model(base_name)

			#output area receptor file
			f_output_receptor_buried_area = os.path.join(path_analysis_pdb_complex_b.value,base_name+".outAreaRecep")
			#ndx files
			#f_ndx = os.path.join(path_analysis_pdb_complex_b.value,base_name+".ndx")
			f_ndx_temporary_index_z = os.path.join(path_analysis_pdb_complex_b.value,base_name+"_temporary_index_z"+".ndx")
			f_ndx_temporary = os.path.join(path_analysis_pdb_complex_b.value,base_name+"_temporary"+".ndx")
			f_ndx_temporary_sasa = os.path.join(path_analysis_pdb_complex_b.value,base_name+"_temporary_sasa"+".ndx")

			#xvg files
			f_xvg_temporary_sasa_res_lig = os.path.join(path_analysis_pdb_complex_b.value,base_name+"_temporary_sasa_res-lig"+".xvg")
			f_xvg_temporary_sasa_res  = os.path.join(path_analysis_pdb_complex_b.value,base_name+"_temporary_sasa_res"+".xvg")
			f_xvg_temporary_sasa_rec_lig  = os.path.join(path_analysis_pdb_complex_b.value,base_name+"_temporary_sasa_rec_lig"+".xvg")
			f_xvg_temporary_sasa_rec  = os.path.join(path_analysis_pdb_complex_b.value,base_name+"_temporary_sasa_rec"+".xvg")

			# Creates a selection with the residues that are closer than 6A to the ligand
			script_make_ndx_buried_area_receptor = SparkFiles.get("make_ndx_buried_area_receptor.sh") #Getting bash script that was copied by addFile command
			command = script_make_ndx_buried_area_receptor + " " + gromacs_path.value + " "+ pdb_complex + " "+ f_ndx_temporary_index_z + " "+ f_ndx_temporary
			process = Popen(command,shell=True, stdout=PIPE, stderr=PIPE)
			stdout, stderr = process.communicate()
			#coping file
			if os.path.exists(f_ndx_temporary):
				shutil.copy(f_ndx_temporary, f_ndx_temporary_sasa)
				#Get all residues for computing area receptor
				all_res = get_residues_receptor_from_ndx_files(f_ndx_temporary)
				returned_list = []
				for res in all_res:
					script_make_ndx_buried_area_receptor_res = SparkFiles.get("make_ndx_buried_area_receptor_res.sh") #Getting bash script that was copied by addFile command
					command = script_make_ndx_buried_area_receptor_res + " " + gromacs_path.value + " "+ pdb_complex + " "+ f_ndx_temporary_sasa + " "+ str(res)
					process = Popen(command,shell=True, stdout=PIPE, stderr=PIPE)
					stdout, stderr = process.communicate()
					# compute surface of system - saved on xvg
					command = gromacs_path.value +"gmx sasa -surface complex -output rec_"+str(res)+ " -o "+ f_xvg_temporary_sasa_res_lig + " -xvg none -f " + pdb_complex +" -s " + pdb_complex + " -n "+ f_ndx_temporary + " -nopbc "
					process = Popen(command,shell=True, stdout=PIPE, stderr=PIPE)
					stdout, stderr = process.communicate()
					# compute surface of receptor - save on xvg
					command = gromacs_path.value +"gmx sasa -surface rec -output rec_"+str(res)+ " -o "+ f_xvg_temporary_sasa_res + " -xvg none -f " + pdb_complex +" -s " + pdb_complex + " -n "+ f_ndx_temporary + " -nopbc "
					process = Popen(command,shell=True, stdout=PIPE, stderr=PIPE)
					stdout, stderr = process.communicate()
					#calculate area
					if os.path.exists(f_xvg_temporary_sasa_res_lig):
						buried_receptor_system = get_value_from_xvg_sasa(f_xvg_temporary_sasa_res_lig)
					else:
						buried_receptor_system = 0
					if os.path.exists(f_xvg_temporary_sasa_res):
						buried_receptor_res  = get_value_from_xvg_sasa(f_xvg_temporary_sasa_res)
					else:
						buried_receptor_res = 0
					res_buried_area = buried_receptor_res - buried_receptor_system
					if (res_buried_area > 0) and (buried_receptor_res > 0):
						res_buried_area_perc = res_buried_area/buried_receptor_res
						#Generating result
						result = (base_name, res, res_buried_area,  res_buried_area_perc)
						returned_list.append(result)
					#Deleting files
					if os.path.exists(f_xvg_temporary_sasa_res_lig):
						os.remove(f_xvg_temporary_sasa_res_lig)
					if os.path.exists(f_xvg_temporary_sasa_res):
						os.remove(f_xvg_temporary_sasa_res)

					#Computing Receptor Area
					command = gromacs_path.value +"gmx sasa -surface complex -output rec"+ " -o "+ f_xvg_temporary_sasa_rec_lig + " -xvg none -f " + pdb_complex +" -s " + pdb_complex + " -n "+ f_ndx_temporary + " -nopbc "
					process = Popen(command,shell=True, stdout=PIPE, stderr=PIPE)
					stdout, stderr = process.communicate()

					command = gromacs_path.value +"gmx sasa -surface rec -output rec"+ " -o "+ f_xvg_temporary_sasa_rec + " -xvg none -f " + pdb_complex +" -s " + pdb_complex + " -n "+ f_ndx_temporary + " -nopbc "
					process = Popen(command,shell=True, stdout=PIPE, stderr=PIPE)
					stdout, stderr = process.communicate()

					if os.path.exists(f_xvg_temporary_sasa_rec_lig):
						sasa_rec_lig = get_value_from_xvg_sasa(f_xvg_temporary_sasa_rec_lig)
					else:
						sasa_rec_lig = 0

					if os.path.exists(f_xvg_temporary_sasa_rec):
						sasa_rec = get_value_from_xvg_sasa(f_xvg_temporary_sasa_rec)
					else:
						sasa_rec = 0

					receptor_area = sasa_rec - sasa_rec_lig

					#Saving result file
					output_receptor_buried_area = open(f_output_receptor_buried_area, "w")
					output_receptor_buried_area.write(str(base_name)+" "+str(receptor_area) +"\n")
					output_receptor_buried_area.close()

					#Deleting all files
					if os.path.exists(f_xvg_temporary_sasa_rec_lig):
						os.remove(f_xvg_temporary_sasa_rec_lig)
					if os.path.exists(f_xvg_temporary_sasa_rec):
						os.remove(f_xvg_temporary_sasa_rec)
					if os.path.exists(f_ndx_temporary):
						os.remove(f_ndx_temporary)
					if os.path.exists(f_ndx_temporary_sasa):
						os.remove(f_ndx_temporary_sasa)
					if os.path.exists(f_ndx_temporary_index_z):
						os.remove(f_ndx_temporary_index_z)

					return returned_list
			else:
				#Here means that some problem for computing area
				return (base_name, "NAN", float(0),  float(0))

# ********** Finish function **********************************************************

# ********** Starting function **********************************************************
		def build_list_model_for_complex(model):
			full_path_model = model[0]
			model_file = model[1]
			path_pdb_complex = path_analysis_pdb_complex_b.value #Obtained from broadcast
			#Building complex file based on model file name
			base_name_model = get_name_model_pdb(full_path_model)
			complex_name = "compl_"+base_name_model+".pdb"
			full_path_for_save_complex = os.path.join(path_pdb_complex,complex_name)
			list_receptor_model_file = (model_file, full_path_for_save_complex)
			save_model_receptor(list_receptor_model_file)
			list_ret = compute_buried_area_all_residues_and_receptor_area(full_path_for_save_complex)
			if os.path.exists(full_path_for_save_complex):
				os.remove(full_path_for_save_complex)
			return list_ret
# ********** Finish function **********************************************************

		#Computing buried area of All-residues and receptor
		all_model_filesRDD = sc.parallelize(all_model_filesRDD)
		all_model_filesRDD = all_model_filesRDD.map(build_list_model_for_complex).collect()
		full_area_file  = os.path.join(path_analysis,base_file_name_receptor+".recepArea")
		save_receptor_buried_area_receptor(full_area_file, all_model_filesRDD)

# ***************** Starting ******************************************/
	#Loading All-residues files
	all_area_file = os.path.join(path_analysis,"*.recepArea")
	buried_areaRDD = sc.textFile(all_area_file).map(loading_lines_from_recepArea_files).collect()
	#Sorting by res_buried_area_perc column
	buried_area_sorted_by_res_buried_area_perc = sorting_buried_area_all_residues(sc, buried_areaRDD)
	buried_area_sorted_by_res_buried_area_perc = buried_area_sorted_by_res_buried_area_perc.map(lambda p: (p.res, p.res_buried_area, p.res_buried_area_perc, p.pose) ).collect() #p.receptor, p.ligand, p.model,

	#Saving buried area file
	path_file_buried_area = os.path.join(path_analysis, "all-residue_buried_areas.dat")
	save_buried_area_recep(path_file_buried_area, buried_area_sorted_by_res_buried_area_perc)

	#Removing all area files
	all_area_files = get_files_recepArea(path_analysis)
	for area_file in all_area_files:
		os.remove(area_file)
# ***************** Finish ******************************************/

# ***************** Starting ******************************************/

	#Loading outAreaRecep files
	all_outAreaRecep_file = os.path.join(path_analysis_pdb,"*.outAreaRecep")
	buried_outAreaRecepRDD = sc.textFile(all_outAreaRecep_file).map(loading_lines_from_outAreaRecep_files).collect()

	buried_outAreaRecepRDD_sort_by_buried_lig_rec = sorting_by_buried_lig_rec(sc, buried_outAreaRecepRDD)
	buried_outAreaRecepRDD_sort_by_buried_lig_rec = buried_outAreaRecepRDD_sort_by_buried_lig_rec.map(lambda p: (p.pose, p.buried_lig_rec) ).collect()

	#Saving buried area receptor file
	path_file_buried_area_rec = os.path.join(path_analysis, "summary_buried_areas_receptor.dat")
	save_buried_area_receptor_sort(path_file_buried_area_rec, buried_outAreaRecepRDD_sort_by_buried_lig_rec)

	#Removing all outAreaRecep files
	all_outAreaRecep = get_files_outAreaRecep(path_analysis)
	for outAreaRecep in all_outAreaRecep:
		os.remove(outAreaRecep)
# ***************** Finish ******************************************/


	finish_time = datetime.now()

	save_log(finish_time, start_time)


main()
