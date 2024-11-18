#!/usr/bin/env python3
#
# Lennard-Jones-Drill-2: CarbonDioxide
# ------------------------------------

# Imports
# -------

import textwrap

class CarbonDioxide(object):

    def __init__(self):

        self.resi_name = 'CO2'

    def get_monomer_a_species(self):

        '''

        Get the Monomer A Species

        '''

        monomer_a_species = {
            'O1': self.get_sp2_oxygen_zmatrix(),
            'C1': self.get_monomer_a_carbon_sp2_zmatrix()
        }

        return monomer_a_species

    def get_sp2_oxygen_zmatrix(self):

        '''

        '''

        zmatrix = '''\
            O11
            C11 O11 1.1802
            O12 C11 1.1802 O11 180.0000
            X11 O11 1.0000 C11 90.0000 O12 180.0000
            0 1
        '''

        atom_name = [
          'O1', 'C1', 'O2'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_a_carbon_sp2_zmatrix(self):

        zmatrix = '''\
            C11
            O11 C11 1.1802
            O12 C11 1.1802 O11 180.0
            X11 C11 1.0000 O11  90.0000 O12 180.0000
            0 1
        '''

        atom_name = [
          'C1', 'O1', 'O2'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_b_species(self):

        monomer_b_species = {
           'O1': self.get_monomer_b_oxygen_zmatrix()
        }

        return monomer_b_species

    def get_monomer_b_oxygen_zmatrix(self):


        zmatrix = '''\
            O21   :1  DISTANCE :2   ANGLE    :3  DIHEDRAL
            X21  O21  1.0000   :1  180.0000   :2   0.0000
            C21 O21 1.1802     :1  180.0000   :2   180.0000
            O22 C21 1.1802    X21 180.0000   O21    0.0000
            0 1
        '''

        atom_name = [
          'O1', 'C1', 'O2'
        ]

        return textwrap.dedent(zmatrix), atom_name
