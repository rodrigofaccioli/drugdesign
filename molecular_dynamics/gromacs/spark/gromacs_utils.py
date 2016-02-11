import os

def check_file_exists(full_path):
	if not os.path.isfile(full_path):
		message = "File NOT found \n"+full_path
		raise Exception(message)