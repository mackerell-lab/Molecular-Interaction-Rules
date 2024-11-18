#!/usr/bin/env python3
#
# Lennard-Jones-Drill-2: Oxetane
# ------------------------------

# Imports
# -------

import textwrap

class Oxetane(object):

    def __init__(self):

        self.resi_name = 'OXTN'

    def get_monomer_a_species(self):

        '''

        Get the Monomer A Species

        '''

        monomer_a_species = {
            'O1': self.get_monomer_a_oxygen_zmatrix(),
        }

        return monomer_a_species

    def get_monomer_a_oxygen_zmatrix(self):

        zmatrix = '''\
          O11
          C11  O11 1.4643
          C12  C11 1.5367 O11  91.6480
          C13  C12 1.5367 C11  83.9135  O11    0.0009
          H11  C13 1.0841 C12 115.7845  C11  115.2619
          H12  C13 1.0841 C12 115.7849  C11 -115.2639
          H13  C12 1.0841 C11 115.4944  O11 -115.3849
          H14  C12 1.0841 C11 115.4945  O11  115.3869
          H15  C11 1.0841 C12 115.7849  C13  115.2639
          H16  C11 1.0841 C12 115.7845  C13 -115.2619
          X11  O11 1.0000 C11  90.0000  C12  180.0000
          0 1
        '''

        atom_name = [
          'O1', 'C1', 'C3', 'C2', 'H7', 'H6', 'H8', 'H9', 'H4', 'H5'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_a_hydrogen_zmatrix(self):

        zmatrix = '''\
          H11
          C11  H11 1.0841
          O11  C11 1.4643 H11 115.7845
          C12  C11 1.5367 O11  91.6480  H11 -115.2619
          C13  C12 1.5367 C11  83.9135  O11    0.0009
          H11  C13 1.0841 C12 115.7845  C11  115.2619
          H12  C13 1.0841 C12 115.7849  C11 -115.2639
          H13  C12 1.0841 C11 115.4944  O11 -115.3849
          H14  C12 1.0841 C11 115.4945  O11  115.3869
          H15  C11 1.0841 C12 115.7849  C13  115.2639
          X11  H11 1.0000 C11  90.0000  O11  180.0000
          0 1
        '''

        atom_name = [
          'H5', 'C1', 'O1', 'C3', 'C2', 'H7', 'H6', 'H8', 'H9', 'H4',
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_b_species(self):

        monomer_b_species = {
            'O1': self.get_monomer_b_oxygen_zmatrix()
        }

        return monomer_b_species

    def get_monomer_b_oxygen_zmatrix(self):

        zmatrix = '''\
            O21  :1  DISTANCE  :2   ANGLE   :3  DIHEDRAL
            C21  O21 1.4643    :1   135.0000   :2    0.0000
            C22  C21 1.5367   O21   90.0000    :1  180.0000
            C23  C22 1.5367 C21  83.9135  O21    0.0000
            H21  C23 1.0841 C22 115.7845  C21  115.2619
            H22  C23 1.0841 C22 115.7849  C21 -115.2639
            H23  C22 1.0841 C21 115.4944  O21 -115.3849
            H24  C22 1.0841 C21 115.4945  O21  115.3869
            H25  C21 1.0841 C22 115.7849  C23  115.2639
            H26  C21 1.0841 C22 115.7845  C23 -115.2619
            X21  O21 1.0000 C21  90.0000  C22  180.0000
            0 1
        '''

        atom_name = [
          'O1', 'C1', 'C3', 'C2', 'H7', 'H6', 'H8', 'H9', 'H4', 'H5'
        ]

        return textwrap.dedent(zmatrix), atom_name

