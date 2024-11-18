#!/usr/bin/env python3
#
# Lennard-Jones-Drill-2: Bromoethane
# ----------------------------------

# Imports
# -------

import textwrap

class Bromoethane(object):

    def __init__(self):

        self.resi_name = 'BRET'

    def get_monomer_a_species(self):

        '''

        Get the Monomer A Species

        '''

        monomer_a_species = {
            'BR1': self.get_monomer_a_bromo_zmatrix()
        }

        return monomer_a_species

    def get_monomer_b_species(self):

      '''

      Get the Monomer A Species

      '''

      monomer_b_species = {
        'BR1': self.get_monomer_b_bromo_zmatrix()
      }

      return monomer_b_species

    def get_monomer_a_bromo_zmatrix(self):

        zmatrix = '''\
            BR11
            C11 BR11 1.9400
            C12 C11 1.5630 BR11 111.6187
            X11 BR11 1.0000 C11  90.0000 C12  180.0000
            H11 C12 1.1042 C11 111.0655 BR11 -59.6837
            H12 C12 1.1042 C11 111.0655 BR11  59.6837
            H13 C12 1.1042 C11 109.5586 BR11 180.0000
            H14 C11 1.1067 C12 110.6939 H11  61.0524
            H15 C11 1.1067 C12 110.6939 H11 180.0000
            0 1
        '''

        atom_name = [
            'BR11', 'C1', 'C2', 'H21', 'H22', 'H23', 'H11', 'H12'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_b_bromo_zmatrix(self):

        zmatrix = '''\
              BR21  :1 DISTANCE  :2 ANGLE  :3  DIHEDRAL
              X21 BR21 1.0000  :1  90.0000   :2   0.0000
              C21 BR21 1.9400 X21  90.0000   :1  180.0000
              C22 C21 1.5630 BR21 111.6187   :1  180.0000
              H21 C22 1.1042 C21 111.0655 BR21 -59.6837
              H22 C22 1.1042 C21 111.0655 BR21  59.6837
              H23 C22 1.1042 C21 109.5586 BR21 180.0000
              H24 C21 1.1067 C22 110.6939 H21  61.0524
              H25 C21 1.1067 C22 110.6939 H21 180.0000
              0 1
          '''

        atom_name = [
          'BR11', 'C1', 'C2', 'H21', 'H22', 'H23', 'H11', 'H12'
        ]

        return textwrap.dedent(zmatrix), atom_name

