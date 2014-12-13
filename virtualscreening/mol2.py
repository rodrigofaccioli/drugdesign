""" 
    Routines for interacting with mol2 file
    These routines were developed by:
    Rodrigo Antonio Faccioli - rodrigo.faccioli@usp.br / rodrigo.faccioli@gmail.com  
    Leandro Oliveira Bortot  - leandro.bortot@usp.br / leandro.obt@gmail.com 
"""

"""
	Returns the molecule name 
	
	pathfilename contains path and file name in which
	wants to know the molecule name.
	In mol2 files, the molecule name is forward line 
	that contains @<TRIPOS>MOLECULE
"""
def get_molecule_name(pathfilename):
	next_line = False
	fmol2 = open(pathfilename, "r")
	for line in fmol2:
		if next_line == True:
			return str(line).strip()
		if line.find("@<TRIPOS>MOLECULE") >= 0:
			next_line = True
