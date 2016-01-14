def adding_chain_atom_line(atom_line,chain="Z"):
	returned_line = ""
	for i in range(0,len(atom_line)):
		#it means column 22
		if i == 21:
			returned_line += chain			
		returned_line += atom_line[i]		
	return returned_line 

def replace_chain_atom_line(atom_line,chain_ref="d",new_chain="Z"):
	return str(atom_line).replace(chain_ref,new_chain)

