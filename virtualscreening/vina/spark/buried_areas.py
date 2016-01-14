from pyspark import SparkContext, SparkConf, SparkFiles
from pyspark.sql import SQLContext, Row
import ConfigParser as configparser
from subprocess import Popen, PIPE
from datetime import datetime
from vina_utils import get_directory_complex_pdb_analysis, get_files_pdb, get_name_model_pdb, get_ligand_from_receptor_ligand_model, get_separator_filename_mode
import os, sys
from os_util import preparing_path
from gromacs_utils import get_value_from_xvg_sasa

def sorting_buried_area(sc, buried_areaRDD):
	sqlCtx = SQLContext(sc)
	buried_areaRDD = sc.parallelize(buried_areaRDD)
	buried_areaRDD = buried_areaRDD.map(lambda p: Row(model=str(p[0]), sasa_lig_min=float(p[1]), sasa_lig_pose=float(p[2]), sasa_lig_complex=float(p[3]), buried_lig_rec_perc=float(p[4]), buried_lig_lig_perc=float(p[5]) ) )
	buried_area_table = sqlCtx.createDataFrame(buried_areaRDD)	
	buried_area_table.registerTempTable("buried_area")

	buried_area_sorted_by_lig_rec_perc = sqlCtx.sql("SELECT * FROM buried_area ORDER BY buried_lig_lig_perc ") 
	return buried_area_sorted_by_lig_rec_perc

def save_buried_area(path_analysis, buried_area_sorted_by_lig_rec_perc):
	path_file_buried_area = os.path.join(path_analysis, "buried_area.txt")
	f_buried_area = open(path_file_buried_area,"w")
	for area in buried_area_sorted_by_lig_rec_perc:
		splited_line = area[0].split("_-_")
		aux_recep = splited_line[0]
		aux_lig = str(splited_line[1])		
		#preparing receptor
		receptor = str(str(aux_recep).replace("u'"," ").replace("compl_", " ")).strip()
		#preparing ligand
		splited_aux_lig = str(aux_lig).split(get_separator_filename_mode())
		ligand = splited_aux_lig[0]
		model = splited_aux_lig[1]
		sasa_lig_min = "{:.4f}".format(area[1])
		sasa_lig_pose = "{:.4f}".format(area[2])
		sasa_lig_complex = "{:.4f}".format(area[3])
		buried_lig_rec_perc = "{:.4f}".format(area[4])
		buried_lig_lig_perc = "{:.4f}".format(area[5])
		line = receptor+"\t"+ligand+"\t"+model+"\t"+str(sasa_lig_min)+"\t"+str(sasa_lig_pose)+"\t"+str(sasa_lig_complex)+"\t"+str(buried_lig_rec_perc)+"\t"+str(buried_lig_lig_perc)+"\n"
		f_buried_area.write(line)	
	f_buried_area.close()

def save_log(finish_time, start_time):
	log_file_name = 'vs_buried_areas.log'
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

	#Adding Python Source file
	#Path for drugdesign project
	path_spark_drugdesign = config.get('DRUGDESIGN', 'path_spark_drugdesign')	
	sc.addPyFile(os.path.join(path_spark_drugdesign,"vina_utils.py"))
	sc.addPyFile(os.path.join(path_spark_drugdesign,"os_util.py"))
	sc.addPyFile(os.path.join(path_spark_drugdesign,"gromacs_utils.py"))

	#Adding bash scripts
	sc.addFile(os.path.join(path_spark_drugdesign,"make_ndx.sh"))

	#Path for Gromacs project
	gromacs_path = preparing_path(config.get('DRUGDESIGN', 'gromacs_path'))
	#Path where PDB ligand are - They are NOT participated in docking
	pdb_ligand_path = config.get('DEFAULT', 'pdb_ligand_path')
	#Path that contains all files for analysis
	path_analysis = config.get('DEFAULT', 'path_analysis')	
	#Path of pdb complex files	
	path_analysis_pdb_complex = get_directory_complex_pdb_analysis(path_analysis)

	#Parameters form command line
	#Indicates probe. Example: 0.14
	probe = float(sys.argv[1])
	#Indicates ndots. Example: 24
	ndots = int(sys.argv[2])

	start_time = datetime.now()

	os.environ["GMX_MAXBACKUP"]="-1"

	#Getting all pdb complex
	all_complex = get_files_pdb(path_analysis_pdb_complex)

	complexRDD = sc.parallelize(all_complex)
	#broadcast
	gromacs_path = sc.broadcast(gromacs_path)
	path_analysis_pdb_complex = sc.broadcast(path_analysis_pdb_complex)
	pdb_ligand_path = sc.broadcast(pdb_ligand_path)
	probe = sc.broadcast(probe)
	ndots = sc.broadcast(ndots)
# ********** Starting function **********************************************************			
	def compute_buried_area(pdb_complex):
		chZ = "chZ"
		buried_lig_rec_perc = -1
		buried_lig_lig_perc = -1
		base_name = get_name_model_pdb(pdb_complex)		
		ligand_name = get_ligand_from_receptor_ligand_model(base_name)
		f_pdb_ligand_no_docking = os.path.join(pdb_ligand_path.value,ligand_name+".pdb")		
		f_ndx = os.path.join(path_analysis_pdb_complex.value,base_name+".ndx")
		f_temp_lig_pose = os.path.join(path_analysis_pdb_complex.value,"lig_pose_"+base_name+".xvg")
		f_temp_lig_complex = os.path.join(path_analysis_pdb_complex.value,"lig_compl_"+base_name+".xvg")
		f_temp_lig_min = os.path.join(path_analysis_pdb_complex.value,"lig_min_"+base_name+".xvg")
		# Makes the index file with the ligand (chain z) and the rest (non chain z)
		script_make_ndx = SparkFiles.get("make_ndx.sh") #Getting bash script that was copied by addFile command
		command = script_make_ndx + " " + gromacs_path.value + " "+ pdb_complex + " "+ f_ndx	
		process = Popen(command,shell=True, stdout=PIPE, stderr=PIPE)
		stdout, stderr = process.communicate()	
		# SASA of the isolated ligand in the pose conformation		
		command = gromacs_path.value +"gmx sasa -f "+pdb_complex+" -s "+pdb_complex+" -n "+f_ndx+" -o "+f_temp_lig_pose+" -surface "+ chZ +" -output " + chZ +" -nopbc -noprot -probe "+ str(probe.value) +" -ndots "+ str(ndots.value)+" -xvg none "
		process = Popen(command,shell=True, stdout=PIPE, stderr=PIPE)
		stdout, stderr = process.communicate()
		# SASA of the complexed ligand in the pose conformation
		command = gromacs_path.value +"gmx sasa -f "+pdb_complex+" -s "+pdb_complex+" -n "+f_ndx+" -o "+f_temp_lig_complex+" -surface "+ "System" +" -output " + chZ +" -nopbc -noprot -probe "+ str(probe.value) +" -ndots "+ str(ndots.value)+" -xvg none "
		process = Popen(command,shell=True, stdout=PIPE, stderr=PIPE)
		stdout, stderr = process.communicate()			
		# SASA of the isolated ligand in its energy-minimized conformation
		command = gromacs_path.value +"gmx sasa -f "+f_pdb_ligand_no_docking+" -s "+f_pdb_ligand_no_docking+" -o "+f_temp_lig_min+" -surface "+ "System" +" -output " + "System" +" -nopbc -noprot -probe "+ str(probe.value) +" -ndots "+ str(ndots.value)+" -xvg none "
		process = Popen(command,shell=True, stdout=PIPE, stderr=PIPE)
		stdout, stderr = process.communicate()
		#Getting values from xvg files
		sasa_lig_min = get_value_from_xvg_sasa(f_temp_lig_min)
		sasa_lig_pose  = get_value_from_xvg_sasa(f_temp_lig_pose)
		sasa_lig_complex  = get_value_from_xvg_sasa(f_temp_lig_complex)
		#Computing buried areas
		buried_lig_rec_perc = (sasa_lig_pose - sasa_lig_complex) / sasa_lig_pose
		buried_lig_lig_perc = (sasa_lig_min - sasa_lig_pose) / sasa_lig_min 
		#Generating result - See column sorting because resultaed file will be created based on this sorting
		returned_list = (base_name, sasa_lig_min, sasa_lig_pose, sasa_lig_complex, buried_lig_rec_perc, buried_lig_lig_perc)
		#Deleting files
		os.remove(f_ndx)
		os.remove(f_temp_lig_complex)
		os.remove(f_temp_lig_pose)
		os.remove(f_temp_lig_min)
		return returned_list
# ********** Finish function **********************************************************					

	#Getting all values
	buried_areaRDD = complexRDD.map(compute_buried_area).collect()
	#Sorting by buried_lig_rec_perc column
	buried_area_sorted_by_lig_rec_perc = sorting_buried_area(sc, buried_areaRDD)
	buried_area_sorted_by_lig_rec_perc = buried_area_sorted_by_lig_rec_perc.map(lambda p: (p.model, p.sasa_lig_min, p.sasa_lig_pose, p.sasa_lig_complex, p.buried_lig_rec_perc, p.buried_lig_lig_perc) ).collect()
	#Saving buried area file
	save_buried_area(path_analysis, buried_area_sorted_by_lig_rec_perc)	

	finish_time = datetime.now()

	save_log(finish_time, start_time)

main()