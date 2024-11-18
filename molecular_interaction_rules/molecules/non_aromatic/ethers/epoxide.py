#!/usr/bin/env python3
#
# Lennard-Jones-Drill-2: Epoxide
# ------------------------------

# Imports
# -------

import textwrap

class Epoxide(object):

    def __init__(self):

        self.resi_name = '1EOX'

    def get_monomer_a_species(self):

        '''

        Get the Monomer A Species

        '''

        monomer_a_species = {
          'O1': self.get_monomer_a_oxygen_zmatrix()
        }

        return monomer_a_species

    def get_monomer_a_oxygen_zmatrix(self):

        zmatrix = '''\
          O11
          C11 O11 1.4535
          C12 O11 1.4535 C11  59.5066
          H11 C12 1.0934 C11 119.3094 O11 -102.7509
          H12 C12 1.0934 C11 119.3094 O11  102.7509
          H13 C11 1.0934 C12 119.3094 H11  154.4981
          H14 C11 1.0934 C12 119.3094 H11    0.0000
          X11 O11 1.0000 C11  90.0000 C12  180.0000
          0 1
        '''

        atom_name = [
            'O1', 'C1', 'C2', 'H3', 'H4', 'H1', 'H2'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_b_species(self):

        monomer_b_species = {
            'O1': self.get_monomer_b_oxygen_zmatrix()
        }

        return monomer_b_species

    def get_monomer_b_oxygen_zmatrix(self):

        zmatrix = '''\
          O21   :1  DISTANCE  :2  ANGLE    :3  DIHEDRAL
          C21  O21  1.4535   :1  150.4257   :2    0.0000
          C22  O21  1.4535   :1  150.4257   :2  180.0000
          H21  C22 1.0934 C21 119.3094 O21 -102.7509
          H22  C22 1.0934 C21 119.3094 O21  102.7509
          H23  C21 1.0934 C22 119.3094 O21  -102.7509
          H24  C21 1.0934 C22 119.3094 O21  102.7509
          0 1
        '''

        atom_name = [
          'O1', 'C1', 'C2', 'H3', 'H4', 'H1', 'H2'
        ]

        return textwrap.dedent(zmatrix), atom_name
