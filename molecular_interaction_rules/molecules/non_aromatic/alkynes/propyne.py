# Imports
# -------

import textwrap

class Propyne(object):

    def __init__(self):

        self.resi_name = 'PRPY'

    def get_monomer_a_species(self):

        monomer_a_species = {
            'C1': self.get_monomer_a_internal_carbon_zmatrix(),
            'H1': self.get_monomer_a_terminal_hydrogen_zmatrix()
        }

        return monomer_a_species


    def get_monomer_b_species(self):

        monomer_b_species = {
            'H1': self.get_monomer_b_terminal_hydrogen_zmatrix()
        }

        return monomer_b_species

    def get_monomer_a_internal_carbon_zmatrix(self):

        zmatrix = '''\
          C12
          C11  C12  1.4721
          C13  C12  1.2335  C11  180.0000
          H11  C13  1.0737  C12  180.0000   C11    8.4872
          H12  C11  1.1004  C12  110.5316   C13  116.9902
          H13  C11  1.1004  C12  110.5316   C13   -3.0097
          H14  C11  1.1004  C12  110.5316   C13 -123.0097
          X11  C12  1.0000  C11   90.0000   C13   0.0000
          0 1
        '''

        atom_name = [
          'C2', 'C1', 'C3', 'H31', 'H11', 'H12', 'H13'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_a_terminal_hydrogen_zmatrix(self):

        zmatrix = '''\
          H11
          C13  H11  1.0737
          C12  C13  1.2335  H11  180.0000
          C11  C12  1.4721  C13  180.0000   H11    8.4872
          H12  C11  1.1004  C12  110.5316   C13  116.9902
          H13  C11  1.1004  C12  110.5316   C13   -3.0097
          H14  C11  1.1004  C13 -110.5316   C12 -123.0097
          X11  H11  1.0000  C13   90.0000   C12  180.0000
          0 1
        '''

        atom_name = [
          'H31', 'C3', 'C2', 'C1', 'H11', 'H12', 'H13'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_b_terminal_hydrogen_zmatrix(self):

        zmatrix = '''\
            H21   :1   DISTANCE   :2  ANGLE   :3  DIHEDRAL
            X21  H21   1.0000     :1  90.0000    :2    0.0000
            C23  H21   1.0737    X21  90.0000    :1  180.0000
            X22  C23   1.0000    H21  90.0000   X21   0.00000
            C22  C23  1.2335  H21  180.0000    :1  180.0000
            C21  C22  1.4721  C23  180.0000   H21    8.4872
            H22  C21  1.1004  C22  110.5316   C23  116.9902
            H23  C21  1.1004  C22  110.5316   C23   3.0097
            H24  C21  1.1004  C23 -110.5316   C22 -123.0097
            0 1
          '''

        atom_name = [
          'H31', 'C3', 'C2', 'C1', 'H11', 'H12', 'H13'
        ]

        return textwrap.dedent(zmatrix), atom_name



