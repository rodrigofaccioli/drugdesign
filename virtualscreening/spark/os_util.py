import os

def preparing_path(ref_path):
	""" 
	This function prepares path for working
	"""		
	last_char = ref_path[len(ref_path)-1]
	if not last_char == os.sep:		
		ref_path += os.sep		
	return ref_path


