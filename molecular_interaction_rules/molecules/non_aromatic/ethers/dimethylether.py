#!/usr/bin/env python3
#
# Lennard-Jones-Drill-2: Dimethylether
# ------------------------------------

# Imports
# -------

import textwrap

class Dimethylether(object):

    def __init__(self):

        self.resi_name = 'DMEE'

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
            C11    O11    1.4238
            H11    C11    1.0978   O11  107.1835
            H12    C11    1.1065   O11  111.0848    H11 -119.4233
            H13    C11    1.1065   O11  111.0847    H11  119.4233
            C12    O11    1.4238   C11  110.6025    H11  179.9995
            H14    C12    1.1065   O11  111.0848    C11  -60.5758
            H15    C12    1.0978   O11  107.1835    C11 -179.9992
            H16    C12    1.1065   O11  111.0847    C11   60.5774
            X11    C11    1.0000   O11 90.0000      C12  0.00000
            0 1
        '''

        atom_name = [
            'O2', 'C1', 'C3', 'H11', 'H12', 'H13', 'H31', 'H32', 'H33'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_b_species(self):

        monomer_b_species = {
            'O1': self.get_monomer_b_oxygen_zmatrix()
        }

        return monomer_b_species

    def get_monomer_b_oxygen_zmatrix(self):

        zmatrix = '''\
            O21   :1  DISTANCE  :2   ANGLE   :3 DIHEDRAL
            C21  O21 1.4643    :1   135.0000   :2    0.0000
            C22  O21 1.4643    :1   135.0000   :2    180.0000
            H21    C21    1.0978   O21  107.1835    C22    0.0000
            H22    C21    1.1065   O21  111.0848    H21 -119.4233
            H23    C21    1.1065   O21  111.0847    H21  119.4233
            H24    C22    1.1065   O21  111.0848    C21  -60.5758
            H25    C22    1.0978   O21  107.1835    C21 -179.9992
            H26    C22    1.1065   O21  111.0847    C21   60.5774
            0 1
        '''

        atom_name = [
          'O2', 'C1', 'C3', 'H11', 'H12', 'H13', 'H31', 'H32', 'H33'
        ]

        return textwrap.dedent(zmatrix), atom_name

