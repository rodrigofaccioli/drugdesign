import json


def create_json_file(output_file_name, data):
    with open(output_file_name, 'w') as outfile:
        json.dump(data, outfile, indent=4, sort_keys=True, separators=(',', ':'))


def create_jsondata_from_docking_output_file(docking_output):
    f_file = open(docking_output, "r")
    data = {}
    for line in f_file:
        splited_line = str(line).split()
        if not len(splited_line) == 0 and splited_line[0].isdigit():
            key = str(splited_line[0])
            value = ', '.join([splited_line[1], splited_line[2], splited_line[3]])
            data[key] = value

    f_file.close()
    return data
