#!/usr/bin/env python3
#
# Lennard-Jones-Drill-2: Furan
# ----------------------------

# Imports
# -------
import textwrap

class Furan(object):

    def __init__(self):

        self.resi_name = 'FURA'

    def get_monomer_a_species(self):

        '''

        Get the Monomer A Species

        '''

        monomer_a_species = {
            'RC1': self.monomer_a_aromatic_zmatrix(),
            'O1': self.monomer_a_oxygen_zmatrix()
        }

        return monomer_a_species

    def get_monomer_b_species(self):

        monomer_b_species = {
            'O1': self.monomer_b_oxygen_zmatrix(),
            'RC1': self.monomer_b_aromatic_zmatrix()
        }

        return monomer_b_species

    def monomer_a_aromatic_zmatrix(self):

        zmatrix = '''\
            X11
            C11  X11  1.2900
            O11  C11  1.3715 X11   60.0000
            C12  O11  1.3723 C11  106.9043 X11    0.0000
            C13  C12  1.3795 O11  110.3901 C11    0.0000
            C14  C12  2.2522 O11   72.5975 C11    0.0000
            H11  C13  1.0882 C12  126.0164 O11 -180.0000
            H12  C14  1.0882 C11  126.0164 O11 -180.0000
            H13  C12  1.0882 O11  115.7876 C11 -180.0000
            H14  C11  1.0882 O11  115.7876 C12 -180.0000
            0 1
        '''

        atom_name = [
            'C1', 'O5', 'C4', 'C3', 'C2', 'H2', 'H3', 'H1', 'H4'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def monomer_a_oxygen_zmatrix(self):

        zmatrix = '''\
            O11
            C11  O11  1.3718
            C12  O11  1.3723 C11  -106.9043
            C13  C12  1.3795 O11  110.3901 C11    0.0000
            C14  C12  2.2522 O11   72.5975 C11    0.0000
            H11  C13  1.0882 C12  126.0164 O11 -180.0000
            H12  C14  1.0882 C11  126.0164 O11 -180.0000
            H13  C12  1.0882 O11  115.7876 C11 -180.0000
            H14  C11  1.0882 O11  115.7876 C12 -180.0000
            X11  O11  1.0000 C11  180.0000 C12    0.0000
            0 1
        '''

        atom_name = [
           'O5','C1', 'C4', 'C3', 'C2', 'H2', 'H3', 'H1', 'H4'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def monomer_b_oxygen_zmatrix(self):

        zmatrix = '''\
          O21   :1   DISTANCE  :2  ANGLE      :3    DIHEDRAL
          C21  O21  1.3723     :1  126.5480   :2     0.0000
          C22  O21  1.3723     :1  126.9043   :2   180.0000
          C23  C21  1.3795   O21  110.3901 C22    0.0000
          C24  C22  1.3795   O21  110.3901  C21    0.0000
          H21  C23  1.0882   C21  130.0157  O21 -180.0000
          H22  C24  1.0882   C22  130.0157  O21 -180.0000
          H23  C22  1.0882   O21  115.8262  C21 -180.0000
          H24  C21  1.0882   O21  115.8262  C22 -180.0000
          X21  O21  1.0000   C21  180.0000  C22    0.0000
          0 1
        '''

        atom_name = [
            'O5', 'C1', 'C4', 'C3', 'C2', 'H2', 'H3', 'H1', 'H4'
        ],

        return textwrap.dedent(zmatrix), atom_name

    def monomer_b_aromatic_zmatrix(self):

      zmatrix = '''\
            X21   :1  DISTANCE  :2   ANGLE   :3   90.0000
            C21  X21  1.2900    :1   90.0000 :2    0.0000
            O21  C21  1.3718 X21   60.0000   :1  DIHEDRAL
            C22  O21  1.3774 C21  106.9043 X21    0.0000
            C23  C22  1.3774 O21  110.4225 C21    0.0000
            C24  C22  2.2522 O21   72.5975 C21    0.0000
            H21  C23  1.0882 C22  126.0157 O21 -180.0000
            H22  C24  1.0882 C21  126.0157 O21 -180.0000
            H23  C22  1.0882 O21  115.8262 C21 -180.0000
            H24  C21  1.0882 O21  115.8262 C22 -180.0000
            0 1
      '''

      atom_name = [
        'C1', 'O5', 'C4', 'C3', 'C2', 'H2', 'H3', 'H1', 'H4'
      ]

      return textwrap.dedent(zmatrix), atom_name
