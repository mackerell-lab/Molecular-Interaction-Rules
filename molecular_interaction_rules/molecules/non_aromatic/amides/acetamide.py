#!/usr/bin/env python3
#
# Lennard-Jones-Drill-2: Acetamide
# --------------------------------

# Imports
# -------

import textwrap

class Acetamide(object):

    def __init__(self):

        self.resi_name = 'ACEM'

    def get_monomer_a_species(self):

        monomer_a_species = {
          'H1': self.get_monomer_a_nitrogen_zmatrix(),
        }

        return monomer_a_species

    def get_monomer_b_species(self):

        monomer_b_species = {
          'H1': self.get_monomer_b_nitrogen_zmatrix()
        }

        return monomer_b_species

    def get_monomer_a_nitrogen_zmatrix(self):

        zmatrix = '''\
            H15
            N11  H15  1.0104
            C11  N11  1.3686   H15  116.5274
            C12  C11  1.5185   N11  116.1617   H15    0.0000
            H11  C12  1.0996   C11  108.5139   N11 -121.4890
            H12  C12  1.0996   C11  108.5134   N11  121.5323
            H13  C12  1.1001   C11  113.2240   N11    0.0216
            O11  C11  1.2378   C12  122.4208   H11   58.5138
            H14  N11  0.9842   C11  120.3664   C12  179.9797
            X11  H15  1.0000   N11   90.0000   C11  180.0000
            0 1
        '''

        atom_name = [
          'HC', 'N', 'C', 'CC', 'HC1', 'HC2', 'HC3', 'O', 'HT',
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_b_nitrogen_zmatrix(self):

        zmatrix = '''\
              H25   :1  DISTANCE   :2  ANGLE  :3  DIHEDRAL
              X21  H25  1.0000     :1   90.0000  :2    0.0000
              N21  H25  1.0104   X21   90.0000   :2  180.0000
              C21  N21  1.3686   H25  116.5274   :1  180.0000
              C22  C21  1.5185   N21  116.1617   H25    0.0000
              H21  C22  1.0996   C21  108.5139   N21 -121.4890
              H22  C22  1.0996   C21  108.5134   N21  121.5323
              H23  C22  1.1001   C21  113.2240   N21    0.0216
              O21  C21  1.2378   C22  122.4208   H21   58.5138
              H24  N21  0.9842   C21  120.3664   C22  179.9797
              X21  H25  1.0000   N21   90.0000   C21  180.0000
              0 1
          '''

        atom_name = [
          'HC', 'N', 'C', 'CC', 'HC1', 'HC2', 'HC3', 'O', 'HT',
        ]

        return textwrap.dedent(zmatrix), atom_name

