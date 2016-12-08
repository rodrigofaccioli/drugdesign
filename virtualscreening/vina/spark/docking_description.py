# Vina docking
class vd_description:

    def __init__(self, receptor, ligand):
        self.receptor = receptor
        self.ligand = ligand

    def get_receptor(self):
        return self.receptor

    def get_ligand(self):
        return self.ligand


# Single docking
class sd_description:

    def __init__(self, pdb_id, ligand):
        self.pdb_id = pdb_id
        self.ligand = ligand

    def get_pdb_id(self):
        return self.pdb_id

    def get_ligand(self):
        return self.ligand