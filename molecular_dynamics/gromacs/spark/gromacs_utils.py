import os

def check_file_exists(full_path):
	if not os.path.isfile(full_path):
		message = "File NOT found \n"+full_path
		raise Exception(message)


def replace_line_in_file(file, search_to, replace_to):
    with open(file, 'r') as input_file:
        for line in input_file:
            if line.strip() == search_to:
                input_file.write(replace_to)
            else:
                input_file.write(line)