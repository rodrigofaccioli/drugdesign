from pyspark import SparkContext, SparkConf, SparkFiles
from pyspark.sql import SQLContext, Row
import ConfigParser as configparser
from subprocess import Popen, PIPE
from datetime import datetime
from vina_utils import get_directory_complex_pdb_analysis, get_files_pdb, get_name_model_pdb, get_ligand_from_receptor_ligand_model, get_separator_filename_mode, get_directory_pdb_analysis, loading_pdb_2_list, get_name_receptor_pdb, get_files_pdb_filter
import os, sys
from os_utils import preparing_path
from gromacs_utils import get_value_from_xvg_sasa
from pdb_io import replace_chain_atom_line
from database_io import load_database

def sorting_buried_area(sc, buried_areaRDD):
	sqlCtx = SQLContext(sc)
	buried_areaRDD = sc.parallelize(buried_areaRDD)
	#buried_areaRDD = buried_areaRDD.map(lambda p: Row(receptor=str(p[0]), ligand=str(p[1]), model=int(p[2]), buried_lig_rec=float(p[3]),  buried_lig_rec_perc=float(p[4]), buried_lig_lig_perc=float(p[5]) ) )
	buried_areaRDD = buried_areaRDD.map(lambda p: Row(pose=str(p[0]), buried_total=float(p[1]) ) )
	buried_area_table = sqlCtx.createDataFrame(buried_areaRDD)
	buried_area_table.registerTempTable("buried_area")

	buried_area_sorted_by_buried_total = sqlCtx.sql("SELECT * FROM buried_area ORDER BY buried_total DESC")  #buried_lig_lig_perc
	return buried_area_sorted_by_buried_total

def save_receptor_buried_area(path_file_buried_area, buried_area_sorted_by_lig_rec_perc):
	f_buried_area = open(path_file_buried_area,"w")
	for area in buried_area_sorted_by_lig_rec_perc:
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
		buried_total = "{:.4f}".format(area[1])
		#line = receptor+"\t"+ligand+"\t"+model+"\t"+str(buried_lig_rec)+"\t"+str(buried_lig_rec_perc)+"\t"+str(buried_lig_lig_perc)+"\n"
		line = pose+"\t"+str(buried_total)+"\n"
		f_buried_area.write(line)
	f_buried_area.close()

def save_buried_area(path_file_buried_area, buried_area_sorted_by_lig_rec_perc):
	f_buried_area = open(path_file_buried_area,"w")
	line = "# buried_area_total[nm2]\tpose"+"\n"
	f_buried_area.write(line)
	for area in buried_area_sorted_by_lig_rec_perc:
		#receptor = area[0]
		#ligand = area[1]
		#model = area[2]
		pose = str(str(area[0]).replace("compl_", " ")).strip()
		buried_total  = "{:.4f}".format(area[1])
		#buried_lig_rec_perc = "{:.4f}".format(area[4])
		#buried_lig_lig_perc = "{:.4f}".format(area[5])
		#line = receptor+"\t"+ligand+"\t"+str(model)+"\t"+str(buried_lig_rec)+"\t"+str(buried_lig_rec_perc)+"\t"+str(buried_lig_lig_perc)+"\n"
		line = str(buried_total)+"\t"+str(pose)+"\n"
		f_buried_area.write(line)
	f_buried_area.close()

def save_normalized_buried_area(path_file_buried_area, full_dataRDD):
	f_buried_area = open(path_file_buried_area,"w")
	line = "# normalized_buried_area_total[nm2]\tpose"+"\n"
	f_buried_area.write(line)
	for area in full_dataRDD.collect():
		pose = str(str(area[0]).replace("compl_", " ")).strip()
		normalized_buried_total  = "{:.4f}".format(area[1])
		line = str(normalized_buried_total)+"\t"+str(pose)+"\n"
		f_buried_area.write(line)
	f_buried_area.close()

def loading_lines_from_area_files(line):
	line_splited = str(line).split()
	#line_ret = ( str(line_splited[0]), str(line_splited[1]), int(line_splited[2]), float(line_splited[3]), float(line_splited[4]), float(line_splited[5]) )
	line_ret = ( str(line_splited[0]), float(line_splited[1]) )
	return line_ret

def get_files_area(mypath):
	only_mol2_file = []
	for root, dirs, files in os.walk(mypath):
		for file in files:
			if file.endswith(".area"):
				f_path = os.path.join(root,file)
				only_mol2_file.append(f_path)
	return only_mol2_file

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

	config = configparser.ConfigParser()
	config.read('config.ini')

	#Path for Gromacs project
	gromacs_path = preparing_path(config.get('DRUGDESIGN', 'gromacs_path'))
	#Path where PDB ligand are - They are NOT participated in docking
	pdb_ligand_path = config.get('DEFAULT', 'pdb_ligand_path')
	#Path that contains all files for analysis
	path_analysis = config.get('DEFAULT', 'path_analysis')
	#Ligand Database file
	ligand_database  = config.get('DEFAULT', 'ligand_database_path_file')
	#Path where all pdb receptor are
	path_receptor_pdb = config.get('DEFAULT', 'pdb_path')
	#Path for saving pdb files of models generated by VS
	path_analysis_pdb = get_directory_pdb_analysis(path_analysis)

	# Create SPARK config
	maxResultSize = str(config.get('SPARK', 'maxResultSize'))
	conf = (SparkConf().set("spark.driver.maxResultSize", maxResultSize))

	# Create context
	sc = SparkContext(conf=conf)
	sqlCtx = SQLContext(sc)

	#Adding Python Source file
	#Path for drugdesign project
	path_spark_drugdesign = config.get('DRUGDESIGN', 'path_spark_drugdesign')
	sc.addPyFile(os.path.join(path_spark_drugdesign,"vina_utils.py"))
	sc.addPyFile(os.path.join(path_spark_drugdesign,"os_utils.py"))
	sc.addPyFile(os.path.join(path_spark_drugdesign,"gromacs_utils.py"))
	sc.addPyFile(os.path.join(path_spark_drugdesign,"pdb_io.py"))
	sc.addPyFile(os.path.join(path_spark_drugdesign,"database_io.py"))
	sc.addPyFile(os.path.join(path_spark_drugdesign,"json_utils.py"))

	#Adding bash scripts
	sc.addFile(os.path.join(path_spark_drugdesign,"make_ndx_buried_area_total.sh"))
	sc.addFile(os.path.join(path_spark_drugdesign,"make_sasa_rec_buried_area_total.sh"))

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

	#Computing Buried areas
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
		def compute_buried_area(pdb_complex):
			chZ = "chZ"

			sasa_complex = -1.0
			sasa_rec = -1.0
			sasa_lig = -1.0
			buried_total = -1.0

			returned_list = []

			try:
				base_name = get_name_model_pdb(pdb_complex)
				ligand_name = get_ligand_from_receptor_ligand_model(base_name)
				f_pdb_ligand_no_docking = os.path.join(pdb_ligand_path.value,ligand_name+".pdb")
				f_ndx = os.path.join(path_analysis_pdb_complex_b.value,base_name+".ndx")

				f_temp_sasa_complex = os.path.join(path_analysis_pdb_complex_b.value,base_name+"_sasa_complex.xvg")
				f_temp_sasa_rec = os.path.join(path_analysis_pdb_complex_b.value,base_name+"_sasa_rec.xvg")
				f_temp_sasa_lig = os.path.join(path_analysis_pdb_complex_b.value,base_name+"_sasa_lig.xvg")

				# Makes the index file with the ligand (chain z) and the rest (non chain z)
				script_make_ndx = SparkFiles.get("make_ndx_buried_area_total.sh") #Getting bash script that was copied by addFile command
				command = script_make_ndx + " " + gromacs_path.value + " "+ pdb_complex + " "+ f_ndx
				process = Popen(command,shell=True, stdout=PIPE, stderr=PIPE)
				stdout, stderr = process.communicate()

				command = gromacs_path.value +"gmx sasa -f " + pdb_complex + " -s " + pdb_complex + " -nopbc " + " -n " + f_ndx + " -surface System " + " -output System "+ " -xvg none " + " -o " + f_temp_sasa_complex
				process = Popen(command,shell=True, stdout=PIPE, stderr=PIPE)
				stdout, stderr = process.communicate()

				# Makes f_temp_sasa_rec file
				script_make_sasa_rec = SparkFiles.get("make_sasa_rec_buried_area_total.sh") #Getting bash script that was copied by addFile command
				command = script_make_sasa_rec + " " + gromacs_path.value + " "+ pdb_complex + " "+ f_ndx + " " + f_temp_sasa_rec
				process = Popen(command,shell=True, stdout=PIPE, stderr=PIPE)
				stdout, stderr = process.communicate()

				command = gromacs_path.value +"gmx sasa -f " + pdb_complex + " -s " + pdb_complex + " -nopbc " + " -n " + f_ndx + " -surface chZ " + " -output chZ "+ " -xvg none " + " -o " +  f_temp_sasa_lig
				process = Popen(command,shell=True, stdout=PIPE, stderr=PIPE)
				stdout, stderr = process.communicate()

				sasa_complex = get_value_from_xvg_sasa(f_temp_sasa_complex)
				sasa_rec = get_value_from_xvg_sasa(f_temp_sasa_rec)
				sasa_lig = get_value_from_xvg_sasa(f_temp_sasa_lig)

				buried_total = sasa_rec + sasa_lig - sasa_complex

				#Generating result - See column sorting because resultaed file will be created based on this sorting
				returned_list = (base_name, buried_total)

			except:
				returned_list = (base_name, float(0))

			#Deleting files
			if os.path.exists(f_ndx):
				os.remove(f_ndx)
			if os.path.exists(f_temp_sasa_complex):
				os.remove(f_temp_sasa_complex)
			if os.path.exists(f_temp_sasa_rec):
				os.remove(f_temp_sasa_rec)
			if os.path.exists(f_temp_sasa_lig):
				os.remove(f_temp_sasa_lig)

			return returned_list
# ********** Finish function **********************************************************

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
			list_ret = compute_buried_area(full_path_for_save_complex)
			os.remove(full_path_for_save_complex)
			return list_ret
# ********** Finish function **********************************************************

		all_model_filesRDD = sc.parallelize(all_model_filesRDD)
		all_model_filesRDD = all_model_filesRDD.map(build_list_model_for_complex).collect()
		#Saving buried area of receptor
		full_area_file  = os.path.join(path_analysis,base_file_name_receptor+".area")
		save_receptor_buried_area(full_area_file, all_model_filesRDD)

	#Loading all area file
	all_area_file = os.path.join(path_analysis,"*.area")
	buried_areaRDD = sc.textFile(all_area_file).map(loading_lines_from_area_files).collect()

	#Sorting by buried_total column
	buried_area_sorted_by_buried_total = sorting_buried_area(sc, buried_areaRDD)
	buried_area_sorted_by_buried_total.cache()
	buried_area_sorted_by_buried_total_LIST = buried_area_sorted_by_buried_total.map(lambda p: (p.pose, p.buried_total) ).collect()

	#Saving buried area file
	path_file_buried_area = os.path.join(path_analysis, "summary_buried_areas_total.dat")
	save_buried_area(path_file_buried_area, buried_area_sorted_by_buried_total_LIST)

	#Calculating normalized buried area
	#Loading database
	rdd_database = load_database(sc, ligand_database)
	#Creating Dataframe
	database_table = sqlCtx.createDataFrame(rdd_database)
	database_table.registerTempTable("database")

	number_pose_ligandRDD = buried_area_sorted_by_buried_total.map(lambda p: Row(buried_total=int(p.buried_total), ligand=get_ligand_from_receptor_ligand_model(p.pose), pose=str(p.pose) ) ).collect()
	number_pose_ligand_table = sqlCtx.createDataFrame(number_pose_ligandRDD)
	number_pose_ligand_table.registerTempTable("buried_area_total_sort")

	sql = """
			SELECT pose, (b.buried_total / a.heavyAtom) as normalized_buried_area
			FROM database a
			JOIN buried_area_total_sort b ON b.ligand = a.ligand
			ORDER BY normalized_buried_area DESC
	      """
	#Getting all data
	full_dataRDD = sqlCtx.sql(sql)

	#Saving normalized buried area file
	path_file_buried_area = os.path.join(path_analysis, "summary_normalized_buried_areas.dat")
	save_normalized_buried_area(path_file_buried_area, full_dataRDD)

	#Removing all area files
	all_area_files = get_files_area(path_analysis)
	for area_file in all_area_files:
		os.remove(area_file)

	finish_time = datetime.now()

	save_log(finish_time, start_time)

main()
