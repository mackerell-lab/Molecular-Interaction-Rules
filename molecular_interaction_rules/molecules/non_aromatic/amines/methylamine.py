#!/usr/bin/env python3
#
# Lennard-Jones-Drill-2: Methylamine
# ----------------------------------

# Imports
# -------

import textwrap

class Methylamine(object):

    def __init__(self):

        self.resi_name = 'MAM1'

    def get_monomer_a_species(self):

        '''

        Get the Monomer A Species

        '''

        monomer_a_species = {
            'H1': self.get_monomer_a_carbon_zmatrix(),
            'H2': self.get_monomer_a_nitrogen_zmatrix()
        }

        return monomer_a_species

    def get_monomer_a_nitrogen_zmatrix(self):

        zmatrix = '''\
            H11
            N11 H11 1.0199
            C11 N11 1.4728 H11 110.0093
            H12 C11 1.1057 N11 115.0533 H11 150.000
            H13 C11 1.0999 N11 108.8647 H11 30.0000
            H14 C11 1.0999 N11 108.8647 H11 270.0000
            H15 N11 1.0199 C11 110.0093 H12 210.0000
            X11 H11 1.0000 N11 90.0000 C11 180.0000
            0 1
        '''

        atom_name = [
          'HN2', 'N1', 'C1', 'HC1', 'HC2', 'HC3', 'HN1',
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_a_carbon_zmatrix(self):

        zmatrix = '''\
            H11
            C11 H11 1.1057
            N11 C11 1.4728 H11 115.0533
            H12 C11 1.0999 N11 108.8647 H11  121.6124
            H13 C11 1.0999 N11 108.8647 H11 -121.6124
            H14 N11 1.0199 C11 110.0093 H11   58.0867
            H15 N11 1.0199 C11 110.0093 H11  -58.0867
            X11 H11 1.0000 C11  90.0000 N11  180.0000
            0 1
        '''

        atom_name = [
          'HC1', 'C1', 'N1', 'HC2', 'HC3', 'HN1', 'HN2'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_b_species(self):

        monomer_b_species = {
            'H1': self.get_monomer_b_nitrogen_zmatrix(),
            'H2': self.get_monomer_b_carbon_zmatrix()
        }

        return monomer_b_species

    def get_monomer_b_nitrogen_zmatrix(self):

        zmatrix = '''\
            H21 :1 DISTANCE  :2  ANGLE  :3  DIHEDRAL
            X21 H21 1.0000 :1 90.0000 :2 180.0000
            N21 H21 1.0199 X21 90.0000 :1 180.0000
            C21 N21 1.4728 H21 110.0093 X21 0.0000
            H22 C21 1.1057 N21 115.0533 H21 -58.0867
            H23 C21 1.0999 N21 108.8647 H21 121.6124
            H24 C21 1.0999 N21 108.8647 H21 -121.6124
            H25 N21 1.0199 C21 110.0093 H21 58.0867
            0 1
        '''

        atom_name = [
          'HN2', 'N1', 'C1', 'HC1', 'HC2', 'HC3', 'HN1',
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_b_carbon_zmatrix(self):

        zmatrix = '''\
            H21 :1 DISTANCE :2  ANGLE   :3 DIHEDRAL
            X21 H21 1.0000   :1 90.0000 :2 180.0000
            C21 H21 1.1057   X21 90.0000 :1  0.0000
            N21 C21 1.4728 H21 115.0533  X21 0.0000
            H22 C21 1.0999 N21 108.8647 H21  121.6124
            H23 C21 1.0999 N21 108.8647 H21 -121.6124
            H24 N21 1.0199 C21 110.0093 H21   58.0867
            H25 N21 1.0199 C21 110.0093 H21  -58.0867
            0 1
        '''

        atom_name = [
          'HC1', 'C1', 'N1', 'HC2', 'HC3', 'HN1', 'HN2'
        ]

        return textwrap.dedent(zmatrix), atom_name

