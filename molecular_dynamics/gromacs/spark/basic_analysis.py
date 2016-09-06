from subprocess import Popen, PIPE


def run_basic_analysis(basic_an_data):
    gromacs_path, nonwater_xtc, nonwater_tpr, simulation_prefix,\
            ana_dir, time_dt = basic_an_data

    # RMSD
    command = ("echo C-alpha C-alpha | "
               + gromacs_path
               + "./gmx rms "
               + "-f "
               + nonwater_xtc
               + " -s "
               + nonwater_tpr
               + " -o "
               + ana_dir
               + "/rmsd_ca_"
               + simulation_prefix
               + ".xvg")
    print(command + "HAHAHAHAHA\n\n\n\n")
    proc = Popen(command, shell=True, stdout=PIPE)
    proc.communicate()

    # Sasa
    command = ("echo Protein | "
               + gromacs_path
               + "./gmx sasa"
               + " -f "
               + nonwater_xtc
               + " -s "
               + nonwater_tpr
               + " -o "
               + ana_dir
               + "/sasa_"
               + simulation_prefix
               + ".xvg "
               + " -dt "
               + str(time_dt))
    print(command + "HAHAHAHAHA\n\n\n\n")
    proc = Popen(command, shell=True, stdout=PIPE)
    proc.communicate()

    # Gyrate
    command = ("echo C-alpha | "
               + gromacs_path
               + "./gmx gyrate "
               + " -f "
               + nonwater_xtc
               + " -s "
               + nonwater_tpr
               + " -o "
               + ana_dir
               + "/rg_ca_"
               + str(simulation_prefix)
               + ".xvg")
    print(command + "HAHAHAHAHA\n\n\n\n")
    proc = Popen(command, shell=True, stdout=PIPE)
    proc.communicate()

    # Number of hydrogen bonds as a function of time
    command = ("echo Protein Protein | "
               + gromacs_path
               + "./gmx hbond"
               + " -f "
               + nonwater_xtc
               + " -s "
               + nonwater_tpr
               + " -num "
               + ana_dir
               + "/hb_p-p_"
               + simulation_prefix
               + ".xvg")
    print(command + "HAHAHAHAHA\n\n\n\n")
    proc = Popen(command, shell=True, stdout=PIPE)
    proc.communicate()


# if __name__ == '__main__':

    # nonwater_xtc = "/home/hakory/Execute/teste_script_traj/ns3pro/ \
    # analysis/md_non-water.1.xtc"
    # reference_tpr = "/home/hakory/Execute/teste_script_traj/ns3pro/md.1.tpr"

    # def basic_analysis(gromacs_path, nonwater_xtc, reference_tpr, prefix)
