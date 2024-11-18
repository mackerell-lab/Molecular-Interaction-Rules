# Imports
# -------

import textwrap

class Neopentane(object):

    def __init__(self):

        self.resi_name = 'neop'

    def get_monomer_a_species(self):

        '''

        Get the Monomer A Species

        '''

        monomer_a_species = {
            'H1': self.get_monomer_a_substituent_carbon(),
        }

        return monomer_a_species


    def get_monomer_b_species(self):

      monomer_b_species = {
        'H1': self.get_monomer_b_substituent_carbon(),
      }

      return monomer_b_species

    def get_monomer_a_substituent_carbon(self):

        zmatrix = '''\
            H11
            C11  H11 1.1031
            C12  C11 1.5365 H11 110.8940
            C13  C12 1.5365 C11 109.4712 H11  -60.0000
            C14  C12 1.5365 C13 109.4713 H11  -60.0000
            C15  C12 1.5365 C13 109.4712 H11 -180.0000
            H12  C13 1.1031 C12 110.8940 C11   60.0000
            H13  C13 1.1031 C12 110.8940 C11 -180.0000
            H14  C13 1.1031 C12 110.8940 C11  -60.0000
            H15  C14 1.1031 C12 110.8940 C11  180.0000
            H16  C14 1.1031 C12 110.8940 C11  -60.0000
            H17  C14 1.1031 C12 110.8940 C11   60.0000
            H18  C15 1.1031 C12 110.8940 C11  -60.0002
            H10  C15 1.1031 C12 110.8940 C11   60.0000
            H20  C15 1.1031 C12 110.8940 C11  180.0000
            H110 C11 1.1031 C12 110.8940 C13   60.0000
            H111 C11 1.1031 C12 110.8940 C13  180.0000
            X11  H11 1.0000 C11  90.0000 C12  180.0000
            0 1
        '''

        atom_name = [
          'H23', 'C2', 'C', 'C1', 'C3', 'C4', 'H11', 'H12', 'H13', 'H31', 'H32', 'H33', 'H41', 'H42', 'H43', 'H21', 'H22',
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_a_carbon_fully_substituent_zmatrix(self):

        zmatrix = '''\
            C12
            C11  C12 1.5365
            H11  C11 1.1031 C12 110.8940
            C13  C12 1.5365 C11 109.4712 H11  -60.0000
            C14  C12 1.5365 C13 109.4713 H11  -60.0000
            C15  C12 1.5365 C13 109.4712 H11 -180.0000
            H12  C13 1.1031 C12 110.8940 C11   60.0000
            H13  C13 1.1031 C12 110.8940 C11 -180.0000
            H14  C13 1.1031 C12 110.8940 C11  -60.0000
            H15  C14 1.1031 C12 110.8940 C11  180.0000
            H16  C14 1.1031 C12 110.8940 C11  -60.0000
            H17  C14 1.1031 C12 110.8940 C11   60.0000
            H18  C15 1.1031 C12 110.8940 C11  -60.0002
            H10  C15 1.1031 C12 110.8940 C11   60.0000
            H20  C15 1.1031 C12 110.8940 C11  180.0000
            H110 C11 1.1031 C12 110.8940 C13   60.0000
            H111 C11 1.1031 C12 110.8940 C13  180.0000
            X11  C12 1.0000 C11 90.0000 C14  180.0000
            0 1
        '''

        atom_name = [
          'C', 'C2','H23','C1', 'C3', 'C4', 'H11', 'H12', 'H13', 'H31', 'H32', 'H33', 'H41', 'H42', 'H43', 'H21', 'H22',
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_b_substituent_carbon(self):

        zmatrix = '''\
              H21  :1  DISTANCE  :2 ANGLE :3 DIHEDRAL
              X21 H21 1.0000  :1  90.0000  :2    0.0000
              C21 H21 1.1013 X21  90.0000  :2  180.0000
              C22 C21 1.5345 H21 111.6112  :1  180.0000
              C23  C22 1.5365 C21 109.4712 H21  -60.0000
              C24  C22 1.5365 C23 109.4713 H21  -60.0000
              C25  C22 1.5365 C23 109.4712 H21 -180.0000
              H22  C23 1.1031 C22 110.8940 C21   60.0000
              H23  C23 1.1031 C22 110.8940 C21 -180.0000
              H24  C23 1.1031 C22 110.8940 C21  -60.0000
              H25  C24 1.1031 C22 110.8940 C21  180.0000
              H26  C24 1.1031 C22 110.8940 C21  -60.0000
              H27  C24 1.1031 C22 110.8940 C21   60.0000
              H28  C25 1.1031 C22 110.8940 C21  -60.0002
              H20  C25 1.1031 C22 110.8940 C21   60.0000
              H20  C25 1.1031 C22 110.8940 C21  180.0000
              H210 C21 1.1031 C22 110.8940 C23   60.0000
              H211 C21 1.1031 C22 110.8940 C23  180.0000
              0 1
          '''

        atom_name = [
          'H23', 'C2', 'C', 'C1', 'C3', 'C4', 'H11', 'H12', 'H13', 'H31', 'H32', 'H33', 'H41', 'H42', 'H43', 'H21', 'H22',
        ]

        return textwrap.dedent(zmatrix), atom_name


