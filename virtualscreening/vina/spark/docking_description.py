# Vina docking
class vd_description:

    def __init__(self, receptor, ligand):
        self.receptor = receptor
        self.ligand = ligand

    def get_receptor(self):
        return self.receptor

    def get_ligand(self):
        return self.ligand
