# Imports
# -------

import textwrap

class Propene(object):

    def __init__(self):

        self.resi_name = 'PRPE'

    def get_monomer_a_species(self):

        '''

        Get the Monomer A Z-Matrices

        '''

        monomer_a_species = {
            'H1': self.get_monomer_a_terminal_hydrogen_zmatrix()
        }

        return monomer_a_species

    def get_monomer_b_species(self):

        '''

        Get the Monomer B Z-Matrices

        '''

        monomer_b_species = {
          'H1': self.get_monomer_b_terminal_hydrogen_zmatrix()
        }

        return monomer_b_species

    def get_monomer_a_terminal_hydrogen_zmatrix(self):

        zmatrix = '''\
            H11
            C11 H11 1.0927
            C12 C11 1.3505 H11 121.1236
            C13 C12 1.5059 C11 124.5302 H11  180.0000
            H12 C13 1.1027 C12 110.8855 C11  120.5942
            H13 C13 1.1027 C12 110.8855 C11 -120.5942
            H14 C13 1.1027 C12 110.8855 C11    0.0000
            H15 C12 1.0973 C13 116.8570 C11  180.0000
            H16 C11 1.0927 C12 121.1236 C13    0.0000
            X11 H11 1.0000 C11  90.0000 C12 180.0000
            0 1
        '''

        atom_name = [
           'H11', 'C1', 'C2', 'C3', 'H31', 'H32', 'H33', 'H21','H12'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_b_terminal_hydrogen_zmatrix(self):

        zmatrix = '''\
              H21   :1    DISTANCE  :2 ANGLE   :3  DIHEDRAL
              X21   H21   1.0000    :1 90.0000   :2   0.0000
              C21 H21 1.0927  X21  90.0000    :2  180.0000
              C22 C21 1.3505 H21 121.1236     :1  180.0000
              C23 C22 1.5059 C21 124.5302 H21  180.0000
              H22 C23 1.1027 C22 110.8855 C21  120.5942
              H23 C23 1.1027 C22 110.8855 C21 -120.5942
              H24 C23 1.1027 C22 110.8855 C21    0.0000
              H25 C22 1.0973 C23 116.8570 C21  180.0000
              H26 C21 1.0927 C22 121.1236 C23    0.0000
              X21 H21 1.0000 C21  90.0000 C22 180.0000
              0 1
          '''

        atom_name = [
          'H11', 'C1', 'C2', 'C3', 'H31', 'H32', 'H33', 'H21','H12'
        ]

        return textwrap.dedent(zmatrix), atom_name

