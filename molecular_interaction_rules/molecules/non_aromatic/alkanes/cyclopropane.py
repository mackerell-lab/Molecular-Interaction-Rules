# Imports
# -------

import textwrap

class Cyclopropane(object):

    def __init__(self):

        self.resi_name = 'C3'

    def get_monomer_a_species(self):

        monomer_a_species = {
            'H1': self.get_monomer_a_hydrogen_zmatrix()
        }

        return monomer_a_species

    def get_monomer_b_species(self):

        monomer_b_species = {
          'H1': self.get_monomer_a_hydrogen_zmatrix()
        }

        return monomer_b_species

    def get_monomer_a_hydrogen_zmatrix(self):

        zmatrix = '''\
          H11
          C11  H11  1.0910
          C12  C11  1.5188  H11  117.6811
          C13  C11  1.5188  C12   60.0000  H11  107.6302
          H12  C13  1.0910  C11  117.6811  C12 -107.6302
          H13  C13  1.0910  C11  117.6811  C12  107.6302
          H14  C12  1.0910  C13  117.6811  H12 -144.7394
          H15  C12  1.0910  C13  117.6811  H12    0.0000
          H16  C11  1.0910  C13  117.6811  H12  144.7394
          X11  H11  1.0000  C11   90.0000  C12  180.0000
          0 1
        '''

        atom_name = [
          'H11', 'C1', 'C2', 'C3', 'H31', 'H32', 'H22', 'H21', 'H12'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_b_hydrogen_zmatrix(self):

        zmatrix = '''\
            H21  :1 DISTANCE  :2 ANGLE :3  DIHEDRAL
            X21 H21 1.0000  :1  90.0000   :2   0.0000
            C21 H21 1.1064 X21  90.0000   :1  180.0000
            C22 C21 1.5369 H21 110.3361   :1  180.0000
            C23  C22 1.5345   C21 111.6112   H21  180.0000
            H22  C23  1.0910  C21  117.6811  C22 -107.6302
            H23  C23  1.0910  C21  117.6811  C22  107.6302
            H24  C22  1.0910  C23  117.6811  H22 -144.7394
            H25  C22  1.0910  C23  117.6811  H22    0.0000
            H26  C21  1.0910  C23  117.6811  H22  144.7394
            0 1
        '''

        atom_name = [
          'H11', 'C1', 'C2', 'C3', 'H31', 'H32', 'H22', 'H21', 'H12'
        ]

        return textwrap.dedent(zmatrix), atom_name

