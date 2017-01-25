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



def sorting_buried_area_ligand(sc, buried_areaRDD):
	sqlCtx = SQLContext(sc)
	buried_areaRDD = sc.parallelize(buried_areaRDD)
	buried_areaRDD = buried_areaRDD.map(lambda p: Row(pose=str(p[0]), buried_lig_rec=float(p[1]), buried_lig_rec_perc=float(p[2]), buried_lig_lig=float(p[3]), buried_lig_lig_perc=float(p[4]) ) ) #receptor=str(p[0]), ligand=str(p[1]), model=int(p[2]),
	buried_area_table = sqlCtx.createDataFrame(buried_areaRDD)	
	buried_area_table.registerTempTable("buried_area_ligand")

	buried_area_sorted_by_lig_rec_perc = sqlCtx.sql("SELECT * FROM buried_area_ligand ORDER BY buried_lig_rec DESC") 
	return buried_area_sorted_by_lig_rec_perc

def save_buried_area_ligand(path_file_buried_area, buried_area_sorted_by_res_buried_area_perc):	
	f_buried_area = open(path_file_buried_area,"w")
	for area in buried_area_sorted_by_res_buried_area_perc:
		#splited_line = area[0].split("_-_")
		#aux_recep = splited_line[0]
		#aux_lig = str(splited_line[1])		
		#preparing receptor
		#receptor = str(str(aux_recep).replace("compl_", " ")).strip()
		#preparing ligand
		#splited_aux_lig = str(aux_lig).split(get_separator_filename_mode())
		#ligand = splited_aux_lig[0]
		#model = splited_aux_lig[1]
		pose = area[0]
		buried_lig_rec = "{:.4f}".format(area[1])
		buried_lig_rec_perc = "{:.4f}".format(area[2])
		buried_lig_lig = "{:.4f}".format(area[3]) 
		buried_lig_lig_perc = "{:.4f}".format(area[4]) 
		#line = receptor+"\t"+ligand+"\t"+model+"\t"+str(buried_lig_rec)+"\t"+str(buried_lig_rec_perc)+"\t"+str(buried_lig_lig)+"\t"+str(buried_lig_lig_perc)+"\n"
		line = str(pose)+"\t"+str(buried_lig_rec)+"\t"+str(buried_lig_rec_perc)+"\t"+str(buried_lig_lig)+"\t"+str(buried_lig_lig_perc)+"\n"
		f_buried_area.write(line)			
	f_buried_area.close()


def save_buried_area_ligand_sort(path_file_buried_area, buried_area):	
	f_buried_area = open(path_file_buried_area,"w")
	line = "# buried_area_lig[nm2]\tburied_area_lig[%]\tburied_area_lig-lig[%]\tpose"+"\n"
	f_buried_area.write(line)	
	for area in buried_area:
		#receptor = area[0]
		#ligand = area[1]
		#model = area[2]
		pose = str(str(area[0]).replace("compl_", " ")).strip()
		buried_lig_rec = "{:.4f}".format(area[1])
		buried_lig_rec_perc = "{:.4f}".format(area[2])
		buried_lig_lig = "{:.4f}".format(area[3])
		buried_lig_lig_perc = "{:.4f}".format(area[4])
		#line = receptor+"\t"+ligand+"\t"+str(model)+"\t"+str(buried_lig_rec)+"\t"+str(buried_lig_rec_perc)+"\t"+str(buried_lig_lig)+"\t"+str(buried_lig_lig_perc)+"\n"
		line = str(buried_lig_rec)+"\t"+str(buried_lig_rec_perc)+"\t"+str(buried_lig_lig_perc)+"\t"+str(pose)+"\n"
		f_buried_area.write(line)			

	f_buried_area.close()

def loading_lines_from_ligandArea_files(line):
	line_splited = str(line).split()
	#line_ret = ( str(line_splited[0]), str(line_splited[1]), int(line_splited[2]), float(line_splited[3]), float(line_splited[4]), float(line_splited[5]), float(line_splited[6]) )
	line_ret = ( str(line_splited[0]), 	float(line_splited[1]), float(line_splited[2]), float(line_splited[3]), float(line_splited[4]) )
	return line_ret

def get_files_ligandArea(mypath):
	only_mol2_file = []
	for root, dirs, files in os.walk(mypath):
		for file in files:
			if file.endswith(".ligandArea"):
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
	log_file_name = 'vs_buried_areas_ligand_receptor.log'
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
	sc.addFile(os.path.join(path_spark_drugdesign,"make_ndx_buried_area_ligand.sh"))	

	#Parameters form command line
	#Indicates probe. Example: 0.14
	probe = float(sys.argv[1])
	#Indicates ndots. Example: 24
	ndots = int(sys.argv[2])

	#Broadcast
	path_analysis_pdb_complex_b = sc.broadcast(path_analysis_pdb)
	gromacs_path = sc.broadcast(gromacs_path)	 
	pdb_ligand_path = sc.broadcast(pdb_ligand_path)
	probe = sc.broadcast(probe)
	ndots = sc.broadcast(ndots)

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
		def compute_buried_area_ligand(pdb_complex):
			chZ = "chZ"
			buried_lig_rec_perc = -1.0
			buried_lig_rec = -1.0
			buried_lig_lig = -1.0
			buried_lig_lig_perc = -1.0
			base_name = get_name_model_pdb(pdb_complex)		
			ligand_name = get_ligand_from_receptor_ligand_model(base_name)
			receptor_name = get_receptor_from_receptor_ligand_model(base_name)
			pose = get_model_from_receptor_ligand_model(base_name)						
			pdb_before_vs = os.path.join(pdb_ligand_path.value,ligand_name+".pdb")			
			#ndx files					
			f_ndx = os.path.join(path_analysis_pdb_complex_b.value,base_name+".ndx")			
			#xvg files
			xvg_temp_sasa_lig_pose = os.path.join(path_analysis_pdb_complex_b.value,base_name+"_sasa_lig_pose"+".xvg")
			xvg_temp_sasa_lig_complex = os.path.join(path_analysis_pdb_complex_b.value,base_name+"_sasa_lig_complex"+".xvg")
			xvg_temp_sasa_lig_min = os.path.join(path_analysis_pdb_complex_b.value,base_name+"_sasa_lig_min"+".xvg")
			# Creates a selection with the residues that are closer than 6A to the ligand
			script_make_ndx_buried_area_ligand = SparkFiles.get("make_ndx_buried_area_ligand.sh") #Getting bash script that was copied by addFile command
			command = script_make_ndx_buried_area_ligand + " " + gromacs_path.value + " "+ pdb_complex + " "+ f_ndx + " "+  xvg_temp_sasa_lig_pose + " "+ str(probe.value)  + " "+ str(ndots.value)  + " "+  xvg_temp_sasa_lig_complex  + " "+ pdb_before_vs  + " "+  xvg_temp_sasa_lig_min
			process = Popen(command,shell=True, stdout=PIPE, stderr=PIPE)
			stdout, stderr = process.communicate()			
			try:
				# SASA of the isolated ligand in the pose conformation			
				sasa_lig_pose = get_value_from_xvg_sasa(xvg_temp_sasa_lig_pose)
				# SASA of the complexed ligand in the pose conformation
				sasa_lig_complex = get_value_from_xvg_sasa(xvg_temp_sasa_lig_complex)
				# SASA of the isolated ligand in its energy-minimized conformation. Only for carbohydrates!
				sasa_lig_min = get_value_from_xvg_sasa(xvg_temp_sasa_lig_min)
				# Area of the ligand which is buried in the receptor
				buried_lig_rec = sasa_lig_pose - sasa_lig_complex
				buried_lig_rec_perc = buried_lig_rec / sasa_lig_pose
				# Area of the ligand in the pose conformation which is buried in itself when compared to the energy-minimized conformation
				buried_lig_lig = sasa_lig_min - sasa_lig_pose
				buried_lig_lig_perc = buried_lig_lig / sasa_lig_min
				returned_list = (base_name, buried_lig_rec, buried_lig_rec_perc, buried_lig_lig, buried_lig_lig_perc)

				#Deleting files
				os.remove(f_ndx)			
				os.remove(xvg_temp_sasa_lig_pose)
				os.remove(xvg_temp_sasa_lig_complex)
				os.remove(xvg_temp_sasa_lig_min)

				return returned_list
			except:
				return (base_name, float(0.0), float(0.0), float(0.0), float(0.0))
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
			list_ret = compute_buried_area_ligand(full_path_for_save_complex)			
			os.remove(full_path_for_save_complex)
			return list_ret
# ********** Finish function **********************************************************	

		all_model_filesRDD = sc.parallelize(all_model_filesRDD)
		all_model_filesRDD = all_model_filesRDD.map(build_list_model_for_complex).collect()	
		#Saving buried area of residue receptor
		full_area_file  = os.path.join(path_analysis,base_file_name_receptor+".ligandArea")
		save_buried_area_ligand(full_area_file, all_model_filesRDD)

	#Loading all area file 
	all_area_file = os.path.join(path_analysis,"*.ligandArea")		
	buried_areaRDD = sc.textFile(all_area_file).map(loading_lines_from_ligandArea_files).collect()	

	#Sorting by buried_lig_lig column
	buried_area_sorted_by_buried_lig_rec = sorting_buried_area_ligand(sc, buried_areaRDD)
	buried_area_sorted_by_buried_lig_rec = buried_area_sorted_by_buried_lig_rec.map(lambda p: (p.pose, p.buried_lig_rec, p.buried_lig_rec_perc, p.buried_lig_lig, p.buried_lig_lig_perc) ).collect() #p.receptor, p.ligand, p.model

	#Saving buried area ligand file
	path_file_buried_area = os.path.join(path_analysis, "summary_buried_area_ligand.dat")
	save_buried_area_ligand_sort(path_file_buried_area, buried_area_sorted_by_buried_lig_rec)	

	#Removing all area files
	all_area_files = get_files_ligandArea(path_analysis)
	for area_file in all_area_files:
		os.remove(area_file)

	finish_time = datetime.now()

	save_log(finish_time, start_time)


main()