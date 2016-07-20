from pyspark.sql import Row	

def load_database(sc, ligand_database):

	ligand_database_file = sc.textFile(ligand_database)

	#Spliting score file by \t
	header = ligand_database_file.first() #extract header	
	rdd_database_split = ligand_database_file.filter(lambda x:x !=header).map(lambda line: line.split("\t"))
	rdd_database = rdd_database_split.map(lambda p: Row(ligand=str(p[0]), torsion=int(p[1]), atom=int(p[2]), heavyAtom=int(p[3]), hb_donors=int(p[4]), hb_acceptors=int(p[5]), hb_donors_acceptors=int(p[6])))
	return rdd_database

