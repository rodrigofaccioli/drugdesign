
import os

def get_files_mol2(mypath):
	only_mol2_file = []
	for root, dirs, files in os.walk(mypath):
		for file in files:
			if file.endswith(".mol2"):
				f_path = os.path.join(root,file)
				only_mol2_file.append(f_path)			
	return only_mol2_file

