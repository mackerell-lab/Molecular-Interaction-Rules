#!/usr/bin/env python3
#
# Lennard-Jones-Drill-2: Tetramethylammonium
# ------------------------------------------

# Imports
# -------

import textwrap

class Tetramethylammonium(object):

    def __init__(self):

        self.resi_name = 'NC4'

    def get_monomer_a_species(self):

        '''

        Get the Monomer A Species

        '''

        monomer_a_species = {
            'H1': self.get_monomer_a_hydrogen_zmatrix()
        }

        return monomer_a_species

    def get_monomer_a_hydrogen_zmatrix(self):

        zmatrix = '''\
            H11
            C11  H11 1.1031
            N11  C11 1.4468 H11 110.8940
            C12  N11 1.4468 C11 109.4712 H11  -60.0000
            C13  N11 1.4468 C12 109.4713 H11  -60.0000
            C14  N11 1.4468 C12 109.4712 H11 -180.0000
            H12  C12 1.1031 N11 110.8940 C11   60.0000
            H13  C12 1.1031 N11 110.8940 C11 -180.0000
            H14  C12 1.1031 N11 110.8940 C11  -60.0000
            H15  C13 1.1031 N11 110.8940 C11  180.0000
            H16  C13 1.1031 N11 110.8940 C11  -60.0000
            H17  C13 1.1031 N11 110.8940 C11   60.0000
            H18  C14 1.1031 N11 110.8940 C11  -60.0002
            H10  C14 1.1031 N11 110.8940 C11   60.0000
            H20  C14 1.1031 N11 110.8940 C11  180.0000
            H110 C11 1.1031 N11 110.8940 C12   60.0000
            H111 C11 1.1031 N11 110.8940 C12  180.0000
            X11  H11 1.0000 C11 90.0000  N11  180.0000
            1 1
        '''

        atom_name = [
            'N', 'C1', 'H1', 'C2', 'C3', 'C4', 'H21', 'H22', 'H23',
            'H31', 'H32',' H33', 'H41', 'H42', 'H43', 'H12', 'H13'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_b_species(self):

        monomer_b_species = {
          'H1': self.get_monomer_b_hydrogen_zmatrix()
        }

        return monomer_b_species

    def get_monomer_b_hydrogen_zmatrix(self):

        zmatrix = '''\
              H21  :1  DISTANCE  :2 90.0000   :3 90.0000
              X21 H21 1.0000     :1  90.0000  :2    0.0000
              C21  H21 1.1031   X21  90.0000  :2  180.0000
              N21  C21 1.4468 H21 110.8940    :1  180.0000
              C22  N21 1.4468 C21 109.4712 H21  -60.0000
              C23  N21 1.4468 C22 109.4713 H21  -60.0000
              C24  N21 1.4468 C22 109.4712 H21 -180.0000
              H22  C22 1.1031 N21 110.8940 C21   60.0000
              H23  C22 1.1031 N21 110.8940 C21 -180.0000
              H24  C22 1.1031 N21 110.8940 C21  -60.0000
              H25  C23 1.1031 N21 110.8940 C21  180.0000
              H26  C23 1.1031 N21 110.8940 C21  -60.0000
              H27  C23 1.1031 N21 110.8940 C21   60.0000
              H28  C24 1.1031 N21 110.8940 C21  -60.0002
              H20  C24 1.1031 N21 110.8940 C21   60.0000
              H30  C24 1.1031 N21 110.8940 C21  180.0000
              H210 C21 1.1031 N21 110.8940 C22   60.0000
              H211 C21 1.1031 N21 110.8940 C22  180.0000
              X11  H21 1.0000 C21 90.0000  N21  180.0000
              1 1
          '''

        atom_name = [
          'N', 'C1', 'H1', 'C2', 'C3', 'C4', 'H21', 'H22', 'H23',
          'H31', 'H32',' H33', 'H41', 'H42', 'H43', 'H12', 'H13'
        ]

        return textwrap.dedent(zmatrix), atom_name

