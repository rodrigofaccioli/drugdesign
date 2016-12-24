from pyspark.sql import Row	

def get_torsion_atom_from_database(sqlCtx, ligand_name):
	sql = ""		
	sql = "select torsion, atom from database where ligand = "	
	sql += str(str("'")+ligand_name+str("'"))		
	df = sqlCtx.sql(sql).rdd.map( lambda p: (p.torsion, p.atom) ).collect() 
	return df[0][0],df[0][1]	

