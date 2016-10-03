import sys
import os
from subprocess import Popen, PIPE
from gromacs_utils import check_file_exists

def create_output_name(ana_dir, name_output, prefix, rep_number):
    return os.path.join(ana_dir, "".join([prefix,
                                          name_output,
                                          rep_number,
                                          ".xvg "]))

def run_basic_analysis(basic_an_data):
    gromacs_path, nonwater_xtc, nonwater_tpr, simulation_prefix,\
            ana_dir, time_dt = basic_an_data

    simulation_prefix = simulation_prefix.split(".")

    prefix = simulation_prefix[0]
    rep_number = simulation_prefix[1]

    output_rmsd = create_output_name(ana_dir, "_rmsd_ca.", prefix, rep_number)

    output_sasa = create_output_name(ana_dir, "_sasa.", prefix, rep_number)

    output_gyrate = create_output_name(ana_dir, "_rg_ca.", prefix, rep_number)

    output_hbond = create_output_name(ana_dir, "_hb_p-p.", prefix, rep_number)

    # RMSD
    command = "".join(["echo C-alpha C-alpha | ",
                       gromacs_path,
                       "./gmx rms ",
                       "-f ",
                       nonwater_xtc,
                       " -s ",
                       nonwater_tpr,
                       " -o ",
                       output_rmsd])
    proc = Popen(command, shell=True, stdout=PIPE)
    proc.communicate()

    # Sasa
    command = "".join(["echo Protein | ",
                       gromacs_path,
                       "./gmx sasa",
                       " -f ",
                       nonwater_xtc,
                       " -s ",
                       nonwater_tpr,
                       " -o ",
                       output_sasa,
                       " -dt ",
                       str(time_dt)])
    proc = Popen(command, shell=True, stdout=PIPE)
    proc.communicate()

    # Gyrate
    command = "".join(["echo C-alpha | ",
                       gromacs_path,
                       "./gmx gyrate ",
                       " -f ",
                       nonwater_xtc,
                       " -s ",
                       nonwater_tpr,
                       " -o ",
                       output_gyrate])
    proc = Popen(command, shell=True, stdout=PIPE)
    proc.communicate()

    # Number of hydrogen bonds as a function of time
    command = "".join(["echo Protein Protein | ",
                       gromacs_path,
                       "./gmx hbond",
                       " -f ",
                       nonwater_xtc,
                       " -s ",
                       nonwater_tpr,
                       " -num ",
                       output_hbond])
    proc = Popen(command, shell=True, stdout=PIPE)
    proc.communicate()


if __name__ == '__main__':

    test_ba_file = sys.argv[1]
    check_file_exists(test_ba_file)
    basic_an_data = tuple(open(test_ba_file, "r"))
    run_basic_analysis(basic_an_data)

