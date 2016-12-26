import json


def create_json_file(output_file_name, data):
    with open(output_file_name, 'w') as outfile:
        json.dump(data, outfile, indent=4, sort_keys=True, separators=(',', ':'))


# Create a dictionary from log file. This dictionary will be used to create json
def create_jsondata_from_docking_output_file(docking_output):
    f_file = open(docking_output, "r")
    list_dict = []

    for line in f_file:
        splited_line = str(line).split()
        if not len(splited_line) == 0 and splited_line[0].isdigit():
            data = dict(mode=int(splited_line[0]),
                        affinity=float(splited_line[1]),
                        dist_from=float(splited_line[2]),
                        best_mode=float(splited_line[3]))
            list_dict.append(data)

    f_file.close()
    return list_dict
