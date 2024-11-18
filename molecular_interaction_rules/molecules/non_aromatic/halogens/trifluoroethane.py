#!/usr/bin/env python3
#
# Lennard-Jones-Drill-2: Trifluoroethane
# --------------------------------------

# Imports
# -------

import textwrap

class Trifluoroethane(object):

    def __init__(self):

        self.resi_name = 'TFET'

    def get_monomer_a_species(self):

        '''

        Get the Monomer A Species

        '''

        monomer_a_species = {
            'F1': self.get_monomer_a_fluoro_zmatrix()
        }

        return monomer_a_species

    def get_monomer_a_fluoro_zmatrix(self):

        zmatrix = '''\
            F11
            C11 F11 1.3624
            C12 C11 1.5028 F11 111.9888
            H11 C12 1.0974 C11 108.9446 F11 180.0000
            F12 C11 1.3624 C12 111.9889 H11 -60.0000
            F13 C11 1.3624 C12 111.9888 H11  60.0000
            H12 C12 1.0974 C11 108.9446 F11 -60.0000
            H13 C12 1.0974 C11 108.9446 F11  60.0000
            X11 F11 1.0000 C11  90.0000 C12 180.0000
            0 1
        '''

        atom_name = [
          'F21', 'C1', 'C2', 'H11', 'F22', 'F23', 'H12', 'H13'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_b_species(self):

        monomer_b_species = {
            'F1': self.get_monomer_b_fluoro_zmatrix()
        }

        return monomer_b_species

    def get_monomer_b_fluoro_zmatrix(self):

        zmatrix = '''\
              F21  :1  DISTANCE :2  ANGLE   :3  DIHEDRAL
              X21 F21 1.0000  :1  90.0000   :2   0.0000
              C21 F21 1.3624  X21  90.0000   :1  180.0000
              C22 C21 1.5028 F21 111.9888    :1  180.0000
              H21 C22 1.0974 C21 108.9446 F21 180.0000
              F22 C21 1.3624 C22 111.9889 H21 -60.0000
              F23 C21 1.3624 C22 111.9888 H21  60.0000
              H22 C22 1.0974 C21 108.9446 F21 -60.0000
              H23 C22 1.0974 C21 108.9446 F21  60.0000
              0 1
          '''

        atom_name = [
          'F21', 'C1', 'C2', 'H11', 'F22', 'F23', 'H12', 'H13'
        ]

        return textwrap.dedent(zmatrix), atom_name
