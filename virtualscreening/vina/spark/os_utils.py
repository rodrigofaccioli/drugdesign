import os


def make_directory(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)


def preparing_path(ref_path):
    """
    This function prepares path for working
    """
    last_char = ref_path[len(ref_path)-1]
    if not last_char == os.sep:
        ref_path += os.sep
    return ref_path


def check_file_exists(full_path):
    if not os.path.isfile(full_path):
        message = "File NOT found \n"+full_path
        raise Exception(message)


def time_execution_log(finish_time, start_time, log_name):
    """
    Generate log file
    """
    log_file_name = log_name
    current_path = os.getcwd()
    path_file = os.path.join(current_path, log_file_name)
    log_file = open(path_file, 'w')

    diff_time = finish_time - start_time
    msg = "".join(['Starting ', str(start_time), '\n'])
    log_file.write(msg)
    msg = "".join(['Finishing ', str(finish_time), '\n'])
    log_file.write(msg)
    msg = "".join(['Time Execution (seconds): ', str(diff_time.total_seconds()), '\n'])
    log_file.write(msg)


try:
    import urllib.request
    import shutil

    def download_file(url, file_name):
        with urllib.request.urlopen(url) as response, open(file_name, 'wb') as out_file:
            shutil.copyfileobj(response, out_file)

except ImportError:
    import urllib

    def download_file(url, file_name):
        urllib.urlretrieve(url, file_name)
