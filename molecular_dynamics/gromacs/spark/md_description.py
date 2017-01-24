class md_description:
    def __init__(self, path, prefix_ref, repetion_number, title_output, total_running):
        self.path = path
        self.prefix_ref = prefix_ref
        self.repetion_number = repetion_number
        self.title_output = title_output
        self.total_running = total_running

    # path where MD files are
    def get_path(self):
        return self.path

    # prefix for reference files: md for md.xtc, md.tpr and md.edr
    def get_prefix_ref(self):
        return self.prefix_ref

    # number of MD repetion
    def get_repetion_number(self):
        return self.repetion_number

    # title for output file name
    def get_title_output(self):
        return self.title_output

    # total running time in ps
    def get_total_running(self):
        return self.total_running

    # MD files must be prefix.rep. Ex: md.1, md.2
    def get_simulation_prefix(self):
        return str(self.get_prefix_ref()) + "." + str(self.get_repetion_number())


class pdb2gmx:
    def __init__(self, dir_number, pdb_input):
        self.dir_number = dir_number
        self.pdb_input = pdb_input

    def get_dir_number(self):
        return self.dir_number

    def get_pdb_input(self):
        return self.pdb_input

    def get_pdb_prefix(self):
        return self.pdb_input.strip('.pdb')


class genbox:
    def __init__(self, dir_number, box_size):
        self.dir_number = dir_number
        self.box_size = box_size

    def get_dir_number(self):
        return self.dir_number

    def get_box_size(self):
        return self.box_size


class add_ions:
    def __init__(self, dir_number, physiological_pattern):
        self.dir_number = dir_number
        self.physiological_pattern = physiological_pattern

    def get_dir_number(self):
        return self.dir_number

    def get_phys_pattern_value(self):
        return self.physiological_pattern


class unrestricted_minimization:
    def __init__(self, dir_number, ion_is):
        self.dir_number = dir_number
        self.ion_is = ion_is

    def get_dir_number(self):
        return self.dir_number

    def get_ion_is(self):
        return self.ion_is


class restricted_minimization:
    def __init__(self, dir_number, min_none):
        self.dir_number = dir_number
        self.min_none = min_none

    def get_dir_number(self):
        return self.dir_number

    def get_min_none(self):
        return self.min_none


class equilibration:
    def __init__(self, dir_number, min_all):
        self.dir_number = dir_number
        self.min_all = min_all

    def get_dir_number(self):
        return self.dir_number

    def get_min_all(self):
        return self.min_all


class molecular_dynamic:
    def __init__(self, dir_number, eq_gro):
        self.dir_number = dir_number
        self.eq_gro = eq_gro

    def get_dir_number(self):
        return self.dir_number

    def get_eq_gro(self):
        return self.eq_gro