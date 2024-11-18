#!/usr/bin/env python3
#
# Lennard-Jones-Drill-2: Acetic Acid
# ----------------------------------

# Imports
# -------

import textwrap

class AceticAcid(object):

    def __init__(self):

        self.resi_name = 'ACEH'

    def get_monomer_a_species(self):

        '''

        Get the Monomer A Species

        '''

        monomer_a_species = {
            'O2': self.get_oxygen_alcohol_zmatrix(),
            'H1': self.get_hydrogen_carbon_zmatrix(),
            'O1': self.get_carbonyl_oxygen_zmatrix(),
        }

        return monomer_a_species

    def get_hydrogen_carbon_zmatrix(self):

        zmatrix = '''\
            H11
            C12 H11 1.1183
            C11 C12 1.5044 H11 108.2226
            O12 C11 1.2500 C12 118.4946 H11   45.0000
            O11 C11 1.2500 O12 116.9653 C12  180.0000
            H12 C12 1.1143 C11 106.2755 O11  101.0000
            H13 C12 1.1151 C11 110.3728 O11  -18.0000
            H14 O11 0.9736 C11 122.0000 C12  180.0000
            X11 H11 1.0000 C12  45.0000 C11   180.0000
            0 1
        '''

        atom_name = [
          'H1', 'C2', 'C1', 'O2',  'O1', 'H2', 'H3', 'HO1',
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_oxygen_alcohol_zmatrix(self):

        zmatrix = '''\
            O11
            C11 O11 1.2500
            C12 C11 1.5044 O11 118.4946
            O12 C11 1.2500 O11 116.9653 C12  180.0000
            H11 C12 1.1143 C11 106.2755 O11  101.0000
            H12 C12 1.1151 C11 110.3728 O11  -18.0000
            H13 C12 1.1183 C11 108.2226 O11 -139.0000
            H14 O11 0.9736 C11 122.0000 C12  180.0000
            X11 O11 1.0000 C11  90.0000 C12    0.0000
            0 1
        '''

        atom_name = [
          'O1', 'C2', 'C1', 'O2', 'H1', 'H2', 'H3', 'HO1'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_carbonyl_oxygen_zmatrix(self):

        zmatrix = '''\
            O12
            C11 O12 1.2500
            C12 C11 1.5044 O12 118.4946
            O11 C11 1.2500 O12 116.9653 C12 180.0000
            H11 C12 1.1143 C11 106.2755 O11  101.0000
            H12 C12 1.1151 C11 110.3728 O11  -18.0000
            H13 C12 1.1183 C11 108.2226 O11 -139.0000
            H14 O11 0.9736 C11 122.0000 C12  180.0000
            X11 O12 1.0000 C11  90.0000 C12    0.0000
            0 1
        '''

        atom_name = [
          'O2', 'C1', 'C2', 'O1', 'H1', 'H2', 'H3', 'HO1',
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_b_species(self):

        '''

        Get the Monomer B Species

        '''

        monomer_b_species = {
            'H1': self.get_hydrogen_carbon_zmatrix()
        }

        return monomer_b_species

    def get_monomer_b_hydrogen_zmatrix(self):

        zmatrix = '''\
            H21  :1  DISTANCE  :2 ANGLE :3 DIHEDRAL
            X21 H21 1.0000  :1  90.0000   :2    0.0000
            C22 H21 1.1183 X21  90.0000   :2  180.0000
            C21 C22 1.5044 H21 108.2226   :1  180.0000
            O22 C21 1.2500 C22 118.4946 H21   45.0000
            O21 C21 1.2500 O22 116.9653 C22  180.0000
            H22 C22 1.1143 C21 106.2755 O21  101.0000
            H23 C22 1.1151 C21 110.3728 O21  -18.0000
            H24 O21 0.9736 C21 122.0000 C22  180.0000
            X21 H21 1.0000 C22  45.0000 C21   180.0000
            0 1
        '''

        atom_name = [
          'H1', 'C2', 'C1', 'O2',  'O1', 'H2', 'H3', 'HO1',
        ]

        return textwrap.dedent(zmatrix), atom_name
