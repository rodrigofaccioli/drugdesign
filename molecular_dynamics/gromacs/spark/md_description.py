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
