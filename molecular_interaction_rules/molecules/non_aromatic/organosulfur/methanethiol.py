#!/usr/bin/env python3
#
# Lennard-Jones-Drill-2: MethaneThiol
# -----------------------------------

# Imports
# -------

import textwrap

class MethaneThiol(object):

    def __init__(self):

        self.resi_name = 'MESH'

    def get_monomer_a_species(self):

        '''

        Get the Monomer A Species

        '''

        monomer_a_species = {
            'H1': self.get_hydrogen_sulphur_zmatrix()
        }

        return monomer_a_species

    def get_hydrogen_sulphur_zmatrix(self):

        zmatrix = '''\
          H11
          S11  H11 1.3496
          C11  S11 1.8329 H11   96.6463
          H12  C11 1.0990 S11  106.1319  H11  -180.0000
          H13  C11 1.0990 S11  111.1621  H11  -62.2020
          H14  C11 1.0990 S11  111.1621  H11   62.2020
          X11  H11 1.0000 S11   90.0000  C11  180.0000
          0 1
        '''

        atom_name = [
          'H4', 'S', 'CM', 'H1', 'H2', 'H3'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_b_species(self):

        monomer_b_species = {
            'H1': self.get_monomer_b_hydrogen_sulphur_zmatrix()
        }

        return monomer_b_species

    def get_monomer_b_hydrogen_sulphur_zmatrix(self):

        zmatrix = '''\
            H21   :1  DISTANCE   :2   ANGLE  :3   DIHEDRAL
            X21  H21  1.0000     :1   90.0000  :2   180.0000
            S21  H21  1.3496     :1   180.0000  :2   0.0000
            C21  S21  1.8329    H21   96.6463  :1  90.0000
            H22  C21 1.0990 S21  106.1319  H21  -180.0000
            H23  C21 1.0990 S21  111.1621  H21  -62.2020
            H24  C21 1.0990 S21  111.1621  H21   62.2020
            X21  H21 1.0000 S21   90.0000  C21  180.0000
            0 1
        '''

        atom_name = [
          'H4', 'S', 'CM', 'H1', 'H2', 'H3'
        ]

        return textwrap.dedent(zmatrix), atom_name
