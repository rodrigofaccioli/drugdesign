import sys
from json_utils import create_json_data, create_json_file


def log_to_json(log_file, path_to_save):
    log_file_name = str(log_file)
    splited = log_file_name.split('/')
    json_name = path_to_save + splited[-1].replace('log', 'json')
    key_name = json_name.split('.')[0]
    json_data = create_json_data(log_file)
    json_final_data = {key_name: json_data}
    create_json_file(json_name, json_final_data)

if __name__ == '__main__':
    log_file = sys.argv[1]
    log_to_json(log_file, '')
