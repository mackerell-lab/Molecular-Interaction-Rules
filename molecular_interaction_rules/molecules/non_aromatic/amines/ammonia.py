#!/usr/bin/env python3
#
# Lennard-Jones-Drill-2: Ammonia
# ------------------------------

# Imports
# -------

import textwrap

class Ammonia(object):

    def __init__(self):

        self.resi_name = 'amm1'

    def get_monomer_a_species(self):

        '''

        Get the Monomer A Species

        '''

        monomer_a_species = {
            'H1': self.get_monomer_a_hydrogen_zmatrix()
        }

        return monomer_a_species

    def get_monomer_b_species(self):

        monomer_b_species = {
            'H1': self.get_monomer_b_hydrogen_zmatrix()
        }

        return monomer_b_species

    def get_monomer_a_nitrogen_zmatrix(self):

        zmatrix = '''\
            X11
            N11  X11  1.0000
            H11  N11  1.0202 X11 106.8000
            H12  N11  1.0202 X11 106.8000 H11 120.0000
            H13  N11  1.0202 X11 106.8000 H11 240.0000
            0 1
        '''

        atom_name = [
          'N1', 'H11', 'H12', 'H13'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_a_hydrogen_zmatrix(self):

      zmatrix = '''\
          X11
          H11  X11  1.0000
          N11  H11  1.0202 X11 136.8000
          H12  N11  1.0202 X11 106.8000 H11 -90.0000
          H13  N11  1.0202 X11 106.8000 H11 -210.0000
          0 1
        '''

      atom_name = [
        'H11', 'N1', 'H12', 'H13'
      ]

      return textwrap.dedent(zmatrix), atom_name

    def get_monomer_b_hydrogen_zmatrix(self):

        zmatrix = '''\
            X21   :1  DISTANCE  :2   ANGLE :3   DIHEDRAL
            H21  X21  1.0000    :1   90.0000  :2  0.0000
            N21  H21  1.0202 X21 136.8000     :1   90.0000
            H22  N21  1.0202 X21 106.8000 H21 -90.0000
            H23  N21  1.0202 X21 106.8000 H21 -210.0000
            0 1
        '''

        atom_name = [
          'H11', 'N1', 'H12', 'H13'
        ]

        return textwrap.dedent(zmatrix), atom_name
