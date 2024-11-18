#!/usr/bin/env python3
#
# Lennard-Jones-Drill-2: Acetate
# ------------------------------

# Imports
# -------

import textwrap

class Acetate(object):

    def __init__(self):

        self.resi_name = 'ACET'

    def get_monomer_a_species(self):

        monomer_a_species = {
            'O1': self.get_oxygen_anion_zmatrix(),
        }

        return monomer_a_species

    def get_monomer_a_carbon_zmatrix(self):

        zmatrix = '''\
            C11
            O11 C11 1.2720
            C12 C11 1.5044 O11 118.4946
            O12 C11 1.2720 O11 116.9653 C12  180.0000
            H11 C12 1.1143 C11 106.2755 O11  101.0000
            H12 C12 1.1151 C11 110.3728 O11  -18.0000
            H13 C12 1.1183 C11 108.2226 O11 -139.0000
            X11 O11 1.0000 C11  90.0000 C12    0.0000
            -1 1
        '''

        atom_name = [
            'C2', 'O2', 'C1', 'O1', 'H1', 'H2', 'H3'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_oxygen_anion_zmatrix(self):

        zmatrix = '''\
            O11
            C11 O11 1.2720
            C12 C11 1.5044 O11 118.4946
            O12 C11 1.2720 O11 116.9653 C12  180.0000
            H11 C12 1.1143 C11 106.2755 O11  101.0000
            H12 C12 1.1151 C11 110.3728 O11  -18.0000
            H13 C12 1.1183 C11 108.2226 O11 -139.0000
            X11 O11 1.0000 C11  90.0000 O12    0.0000
            -1 1
        '''

        atom_name = [
          'O2', 'C2', 'C1', 'O1', 'H1', 'H2', 'H3'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_b_species(self):

        monomer_b_species = {
            'O1': self.get_monomer_b_oxygen_zmatrix(),
        }

        return monomer_b_species

    def get_monomer_b_oxygen_zmatrix(self):

        zmatrix = '''\
            O21  :1 DISTANCE  :2  ANGLE   :3  DIHEDRAL
            X21 O21 1.0000     :1  90.0000  :2   0.0000
            C21 O21 1.2720    X21  90.0000   :1  180.0000
            C22 C21 1.5044    O21 118.4946  X21    0.0000
            O22 C21 1.2720 O21 116.9653 C22  180.0000
            H21 C22 1.1143 C21 106.2755 O21  101.0000
            H22 C22 1.1151 C21 110.3728 O21  -18.0000
            H23 C22 1.1183 C21 108.2226 O21 -139.0000
            X21 O21 1.0000 C21  90.0000 O22    0.0000
            -1 1
        '''

        atom_name = [
            'O2', 'C2', 'C1', 'O1', 'H1', 'H2', 'H3'
        ]

        return textwrap.dedent(zmatrix), atom_name

