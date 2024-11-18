# Imports
# -------

import textwrap

class OneThreeDibutene(object):

    def __init__(self):

        self.resi_name = '13DB'

    def get_monomer_a_species(self):

        '''

        Get the Monomer A Z-Matrices

        '''

        monomer_a_species = {
          'H1': self.get_sp2_terminal_hydrogen()
        }

        return monomer_a_species

    def get_monomer_b_species(self):

        '''

        Get the Monomer A Z-Matrices

        '''

        monomer_b_species = {
          'H1': self.get_sp2_monomer_b_terminal_hydrogen()
        }

        return monomer_b_species

    def get_sp2_terminal_hydrogen(self):

        zmatrix = '''\
            H11
            C11 H11 1.0923
            C12 C11 1.3566 H11 121.3761
            C13 C12 1.4644 C11 123.4738 H11 -180.0000
            C14 C13 1.3566 C12 123.4738 C11 -180.0000
            H12 C14 1.0943 C13 120.9363 C12    0.0000
            H13 C14 1.0923 C13 121.3761 C12 -180.0000
            H14 C13 1.0972 C12 117.0408 C11    0.0000
            H15 C12 1.0972 C13 117.0408 C14    0.0000
            H16 C11 1.0943 C12 120.9363 C13    0.0000
            X11 H11 1.0000 C11  90.0000 C12    0.0000
            0 1
        '''

        atom_name = [
          'H11', 'C1', 'C2', 'C3', 'C4', 'H41', 'H42', 'H31', 'H22', 'H12',
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_sp2_monomer_b_terminal_hydrogen(self):

        zmatrix = '''\
              H21    :1    DISTANCE     :2 ANGLE     :3 DIHEDRAL
              X21   H21   1.0000        :1  90.0000    :2  0.0000
              C21   H21  1.1066  X21  90.0000      :2  180.0000
              C22   C21  1.5353  H21  110.8679     :1  180.0000
              C23 C22 1.4644 C21 123.4738 H21 -180.0000
              C24 C23 1.3566 C22 123.4738 C21 -180.0000
              H22 C24 1.0943 C23 120.9363 C22    0.0000
              H23 C24 1.0923 C23 121.3761 C22 -180.0000
              H24 C23 1.0972 C22 117.0408 C21    0.0000
              H25 C22 1.0972 C23 117.0408 C24    0.0000
              H26 C21 1.0943 C22 120.9363 C23    0.0000
              0 1
          '''

        atom_name = [
          'H11', 'C1', 'C2', 'C3', 'C4', 'H41', 'H42', 'H31', 'H22', 'H12',
        ]

        return textwrap.dedent(zmatrix), atom_name


