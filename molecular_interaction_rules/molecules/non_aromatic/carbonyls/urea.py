#!/usr/bin/env python3
#
# Lennard-Jones-Drill-2: Urea
# ---------------------------

# Imports
# -------

import textwrap

class Urea(object):

    def __init__(self):

        self.resi_name = ''

    def get_monomer_a_species(self):

        monomer_a_species = {
            'N1': self.get_nitrogen_zmatrix(),
            'O1': self.get_oxygen_zmatrix()
        }

        return monomer_a_species

    def get_nitrogen_zmatrix(self):

        zmatrix = '''\
            H14
            N11 H14 1.0143
            C11 N11 1.3945 H14 112.4245
            N12 C11 1.3945 N11 113.4670 H14    0.0000
            H11 N12 1.0143 C11 112.4245 N11  180.0000
            H12 N12 1.0143 C11 112.4245 N11    0.0000
            O11 C11 1.2296 N12 123.2665 H11    0.0000
            H13 N11 1.0143 C11 112.4245 N12  180.0000
            X11 H14 1.0000 N11  90.0000 C11  180.0000
            0 1
        '''

        atom_name = [
            'H32', 'N1', 'C2', 'N3', 'H11', 'H12', 'O1', 'H31',
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_oxygen_zmatrix(self):

        zmatrix = '''\
            O11
            C11 O11 1.2296
            N11 C11 1.3945 O11  123.2665
            N12 C11 1.3945 N11  113.4670 O11 180.0000
            H11 N12 1.0143 C11 112.4245 N11  180.0000
            H12 N12 1.0143 C11 112.4245 N11    0.0000
            H13 N11 1.0143 C11 112.4245 N12  180.0000
            H14 N11 1.0143 C11 112.4245 N12    0.0000
            X11 O11 1.0000 C11  90.0000 N11  180.0000
            0 1
        '''

        atom_name = [
            'O1', 'C2', 'N1', 'N3', 'H11', 'H12', 'H31', 'H32'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_b_species(self):

        monomer_b_species = {
            'N1': self.get_monomer_b_nitrogen_zmatrix(),
            'O1': self.get_monomer_b_oxygen_zmatrix()
        }

        return monomer_b_species

    def get_monomer_b_nitrogen_zmatrix(self):

        zmatrix = '''\
            H24  :1  DISTANCE   :2  ANGLE     :3  90.0000
            X21  H24  1.00000   :1  90.0000   :2  180.0000
            N21  H24  1.1043   X21 90.0000    :1   180.0000
            C21 N21 1.3945 H24 112.4245  X21 DIHDERAL
            N22 C21 1.3945 N21 113.4670  X21     0.0000
            H21 N22 1.0143 C21 112.4245 N21  180.0000
            H22 N22 1.0143 C21 112.4245 N21    0.0000
            O21 C21 1.2296 N22 123.2665 H21    0.0000
            H23 N21 1.0143 C21 112.4245 N22  180.0000
            H24 N21 1.0143 C21 112.4245 N22    0.0000
            X21 N21 1.0000 C21  90.0000 N22  180.0000
            0 1
        '''

        atom_name = [
           'N1', 'C2', 'N3', 'H11', 'H12', 'O1', 'H31', 'H32'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_b_oxygen_zmatrix(self):

        zmatrix = '''\
            O21  :1  DISTANCE   :2  ANGLE  :3    DIHEDRAL
            C21 O21 1.2296      :1  180.0000 :2   180.0000
            N21 C21 1.3945 O21  123.2665     :1     0.0000
            N22 C21 1.3945 N21  113.4670 O21 180.0000
            H21 N22 1.0143 C21 112.4245 N21  180.0000
            H22 N22 1.0143 C21 112.4245 N21    0.0000
            H23 N21 1.0143 C21 112.4245 N22  180.0000
            H24 N21 1.0143 C21 112.4245 N22    0.0000
            X21 O21 1.0000 C21  90.0000 N21  180.0000
            0 1
        '''

        atom_name = [
          'O1', 'C2', 'N1', 'N3', 'H11', 'H12', 'H31', 'H32'
        ]

        return textwrap.dedent(zmatrix), atom_name

