#!/usr/bin/env python3
#
# Lennard-Jones-Drill-2: DimethylSulfone
# --------------------------------------

# Imports
# -------

import textwrap

class DimethylSulfone(object):

    def __init__(self):

        self.resi_name = 'DMSN'

    def get_monomer_a_species(self):

        '''

        Get the Monomer A Species

        '''

        monomer_a_species = {
            'O1': self.get_oxygen_zmatrix(),
            'S1': self.get_sulphur_matrix()
        }

        return monomer_a_species

    def get_oxygen_zmatrix(self):

        zmatrix = '''\
          O11
          S11 O11 1.4948
          C11 S11 1.7995 O11  107.8073
          H11 C11 1.0984 S11  108.6706 O11 -175.3736
          H12 C11 1.0984 S11  108.6706 O11  -52.7697
          H13 C11 1.0984 S11  105.2498 O11   65.9283
          C12 S11 1.7995 C11  103.5586 H11  -61.3019
          H14 C12 1.0984 S11  108.6706 O11   52.7696
          H15 C12 1.0984 S11  105.2498 O11  -65.9283
          H16 C12 1.0984 S11  108.6706 O11  175.3736
          O12 S11 1.4948 C11  107.8074 H11   52.7696
          X11 O11 1.0000 S11   90.0000 C11  180.0000
          0 1
        '''

        atom_name = [
          'O1', 'S', 'C3', 'H31', 'H32', 'H33', 'C4', 'H41', 'H42', 'H43', 'O2'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_sulphur_matrix(self):

        zmatrix = '''\
          S11
          O11 S11 1.4948
          C11 S11 1.7995 O11  107.8073
          H11 C11 1.0984 S11  108.6706 O11 -175.3736
          H12 C11 1.0984 S11  108.6706 O11  -52.7697
          H13 C11 1.0984 S11  105.2498 O11   65.9283
          C12 S11 1.7995 C11  103.5586 H11  -61.3019
          H14 C12 1.0984 S11  108.6706 O11   52.7696
          H15 C12 1.0984 S11  105.2498 O11  -65.9283
          H16 C12 1.0984 S11  108.6706 O11  175.3736
          O12 S11 1.4948 C11  107.8074 H11   52.7696
          X11 S11 1.0000 O11  55.0000 C11  135.0000
          0 1
        '''

        atom_name = [
          'S', 'O1','C3', 'H31', 'H32', 'H33', 'C4', 'H41', 'H42', 'H43', 'O2'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_b_species(self):

        '''

        Get the Monomer B Species

        '''

        monomer_b_species = {
            'O1': self.get_monomer_b_oxygen_zmatrix(),
            'S1': self.get_monomer_b_sulphur_zmatrix()
        }

        return monomer_b_species

    def get_monomer_b_oxygen_zmatrix(self):

        zmatrix = '''\
            O21  :1  DISTANCE   :2  ANGLE    :3    DIHEDRAL
            S21 O21 1.4948      :1  180.0000  :2   180.0000
            C21 S21 1.7995 O21  107.8073      :1     0.0000
            H21 C21 1.0984 S21  108.6706 O21 -175.3736
            H22 C21 1.0984 S21  108.6706 O21  -52.7697
            H23 C21 1.0984 S21  105.2498 O21   65.9283
            C22 S21 1.7995 C21  103.5586 H21  -61.3019
            H24 C22 1.0984 S21  108.6706 O21   52.7696
            H25 C22 1.0984 S21  105.2498 O21  -65.9283
            H26 C22 1.0984 S21  108.6706 O21  175.3736
            O22 S21 1.4948 C21  107.8074 H21   52.7696
            X21 O21 1.0000 S21   90.0000 C21  180.0000
            0 1
          '''

        atom_name = [
          'O1', 'S', 'C3', 'H31', 'H32', 'H33', 'C4', 'H41', 'H42', 'H43', 'O2'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_b_sulphur_zmatrix(self):

        zmatrix = '''\
              S21  :1  DISTANCE   :2  ANGLE    :3    DIHEDRAL
              O21 S21 1.4948      :1  180.0000  :2   180.0000
              C21 S21 1.7995 O21  107.8073      :1     0.0000
              H21 C21 1.0984 S21  108.6706 O21 -175.3736
              H22 C21 1.0984 S21  108.6706 O21  -52.7697
              H23 C21 1.0984 S21  105.2498 O21   65.9283
              C22 S21 1.7995 C21  103.5586 H21  -61.3019
              H24 C22 1.0984 S21  108.6706 O21   52.7696
              H25 C22 1.0984 S21  105.2498 O21  -65.9283
              H26 C22 1.0984 S21  108.6706 O21  175.3736
              O22 S21 1.4948 C21  107.8074 H21   52.7696
              X21 O21 1.0000 S21   90.0000 C21  180.0000
              0 1
            '''

        atom_name = [
          'S', 'O1', 'C3', 'H31', 'H32', 'H33', 'C4', 'H41', 'H42', 'H43', 'O2'
        ]

        return textwrap.dedent(zmatrix), atom_name
