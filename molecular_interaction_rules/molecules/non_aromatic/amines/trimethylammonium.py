#!/usr/bin/env python3
#
# Lennard-Jones-Drill-2: Trimethylammonium
# ----------------------------------------

# Imports
# -------

import textwrap

class Trimethylammonium(object):

    def __init__(self):

        self.resi_name = 'NC3'

    def get_monomer_a_species(self):

        '''

        Get the Monomer A Species

        '''

        monomer_a_species = {
            'H1': self.get_monomer_a_nitrogen_hydrogen_zmatrix()
        }

        return monomer_a_species

    def get_monomer_a_nitrogen_hydrogen_zmatrix(self):

        zmatrix = '''\
            H20
            N11  H20  1.0841
            C11  N11  1.4451  H20  109.8360
            H11  C11  1.0841  N11  109.8360  H20   90.4316
            H12  C11  1.0841  N11  109.8360  H11 -118.5543
            H13  C11  1.0949  N11  113.0081  H11  120.7228
            C12  N11  1.4451  C11  111.9391  H11  -57.4136
            H14  C12  1.0841  N11  109.8361  C11   57.4137
            H15  C12  1.0841  N11  109.8360  C11  175.9680
            H16  C12  1.0949  N11  113.0081  C11  -63.3091
            C13  N11  1.4451  C11  111.9391  H11  175.9680
            H17  C13  1.0841  N11  109.8360  C11 -175.9681
            H18  C13  1.0841  N11  109.8360  C11  -57.4137
            H19  C13  1.0949  N11  113.0081  C11   63.3090
            X11  H20  1.0000  N11   90.0000  C11    0.0000
            1 1
        '''

        atom_name = [
          'HN4', 'N', 'C1', 'H11', 'H12', 'H13', 'C2', 'H21', 'H22', 'H23', 'C3', 'C31', 'C32', 'C33'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_b_species(self):

        monomer_b_species = {
            'H1': self.get_monomer_b_nitrogen_hydrogen_zmatrix()
        }

        return monomer_b_species

    def get_monomer_b_nitrogen_hydrogen_zmatrix(self):

        zmatrix = '''\
              H30  :1  DISTANCE   :2   180.0000  :3   180.0000
              N21  H30  1.0841    :1   0.0000    :2   180.0000
              C21  N21  1.4451  H30  109.8360    :1    0.0000
              H21  C21  1.0841  N21  109.8360  H30   90.4316
              H22  C21  1.0841  N21  109.8360  H21 -118.5543
              H23  C21  1.0949  N21  113.0081  H21  120.7228
              C22  N21  1.4451  C21  111.9391  H21  -57.4136
              H24  C22  1.0841  N21  109.8361  C21   57.4137
              H25  C22  1.0841  N21  109.8360  C21  175.9680
              H26  C22  1.0949  N21  113.0081  C21  -63.3091
              C23  N21  1.4451  C21  111.9391  H21  175.9680
              H27  C23  1.0841  N21  109.8360  C21 -175.9681
              H28  C23  1.0841  N21  109.8360  C21  -57.4137
              H29  C23  1.0949  N21  113.0081  C21   63.3090
              X21  H30  1.0000  N21   90.0000  C21    0.0000
              1 1
          '''

        atom_name = [
          'HN4', 'N', 'C1', 'H11', 'H12', 'H13', 'C2', 'H21', 'H22', 'H23', 'C3', 'C31', 'C32', 'C33'
        ]

        return textwrap.dedent(zmatrix), atom_name

