#!/usr/bin/env python3
#
# Lennard-Jones-Drill-2: MethylAcetate
# ------------------------------------

# Imports
# -------

import textwrap

class MethylAcetate(object):

    def __init__(self):

        self.resi_name = 'MAS'

    def get_monomer_a_species(self):

        monomer_a_species = {
            'O1': self.get_monomer_a_oxygen_zmatrix()
        }

        return monomer_a_species

    def get_monomer_a_oxygen_zmatrix(self):

        zmatrix = '''\
            O11
            C11   O11  1.3264
            C12   C11  1.5042   O11  111.4090
            H11   C12  1.0796   C11  109.5903    O11  179.9995
            H12   C12  1.0838   C11  109.7272    O11  -59.0373
            H13   C12  1.0838   C11  109.7273    O11   59.0363
            O11   C11  1.1880   C12  125.1960    H11   -0.0004
            C13   O11  1.4164   C11  116.9354    C12  179.9999
            H14   C13  1.0804   O11  110.5624    C11  -60.5387
            H15   C13  1.0787   O11  105.8609    C11 -179.9995
            H16   C13  1.0804   O11  110.5624    C11   60.5397
            X11   O11  1.0000   C11   90.0000    C12    0.0000
            0 1
        '''

        atom_name = [
          'OM', 'C', 'C2', 'H21', 'H22', 'H23', 'O', 'C1', 'H11', 'H12', 'H13'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_b_species(self):

      '''

      Get the Monomer B Z-Matrices

      '''

      monomer_b_species = {
          'O1': self.get_monomer_b_oxygen_zmatrix()
      }

      return monomer_b_species

    def get_monomer_b_oxygen_zmatrix(self):

      zmatrix = '''\
              O21    :1  DISTANCE  :2  ANGLE      :3    DIHEDRAL
              C21   O21  1.3264    :1  126.2746   :2     0.0000
              C22   C21  1.5042   O21  126.2746   :2    180.0000
              H21   C22  1.0796   C21  109.5903    O21  179.9995
              H22   C22  1.0838   C21  109.7272    O21  -59.0373
              H23   C22  1.0838   C21  109.7273    O21   59.0363
              O21   C21  1.1880   C22  125.1960    H21   -0.0004
              C23   O21  1.4164   C21  116.9354    C22  179.9999
              H24   C23  1.0804   O21  110.5624    C21  -60.5387
              H25   C23  1.0787   O21  105.8609    C21 -179.9995
              H26   C23  1.0804   O21  110.5624    C21   60.5397
              X21   O21  1.0000   C21   90.0000    C22    0.0000
              0 1
      '''

      atom_name = [
          'OM', 'C', 'C2', 'H21', 'H22', 'H23', 'O', 'C1', 'H11', 'H12', 'H13'
      ]

      return textwrap.dedent(zmatrix), atom_name
