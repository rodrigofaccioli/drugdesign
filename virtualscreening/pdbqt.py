""" 
    Routines for manipulating pdbqt files
    These routines were developed by:
    Rodrigo Antonio Faccioli - rodrigo.faccioli@usp.br / rodrigo.faccioli@gmail.com  
    Leandro Oliveira Bortot  - leandro.bortot@usp.br / leandro.obt@gmail.com 
"""

import operator

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
	return number

def sort_dictionary(dic_docking):
	return sorted(dic_docking.items(), key=operator.itemgetter(1), reverse=True)
