#!/usr/bin/env python3
#
# Lennard-Jones-Drill-2: Acetonitrile
# -----------------------------------

# Imports
# -------

import textwrap

class Acetonitrile(object):

    def __init__(self):

        self.resi_name = 'ACN'

    def get_monomer_a_species(self):

        '''

        Get the Monomer A Species

        '''

        monomer_a_species = {
            'N1': self.monomer_a_nitrogen_z_matrix(),
            'C1': self.monomer_a_carbon_z_matrix()
        }

        return monomer_a_species

    def get_monomer_b_species(self):

        monomer_b_species = {
            'N1': self.monomer_b_nitrogen_z_matrix(),
            'C1': self.monomer_b_carbon_z_matrix()
        }

        return monomer_b_species

    def monomer_a_carbon_z_matrix(self):

        zmatrix = '''\
            C11
            N11 C11 1.1854
            X11 N11 1.0000 C11  90.0000
            C12 C11 1.4710 X11  90.0000 C11  0.0000
            H11 C12 1.0990 C11 109.8381 N11 279.5961
            H12 C12 1.0990 C11 109.8381 N11  20.4038
            H13 C12 1.0990 C11 109.8381 N11 140.4038
            X11 C11 1.0000 N11  90.0000 C12 180.0000
            0 1
        '''

        atom_name = [
          'C2', 'N3','C1', 'H11', 'H12', 'H13'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def monomer_a_nitrogen_z_matrix(self):

        zmatrix = '''\
            N11
            C11 N11 1.1854
            C12 C11 1.4710 N11 180.0000
            H11 C12 1.0990 C11 109.8381 N11 279.5961
            H12 C12 1.0990 C11 109.8381 N11  20.4038
            H13 C12 1.0990 C11 109.8381 N11 140.4038
            X11 N11 1.0000 C11  90.0000 C12 180.0000
            0 1
        '''

        atom_name = [
            'N3', 'C2', 'C1', 'H11', 'H12', 'H13'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def monomer_b_nitrogen_z_matrix(self):

        zmatrix = '''\
              N21  :1   DISTANCE  :2  ANGLE    :3  DIHEDRAL
              X21 N21   1.0000    :1  90.0000    :2    0.0000
              C21 N21   1.1854   X21  90.0000    :1  180.0000
              X22 C21   1.0000   N21  90.0000   X21   0.00000
              C22 C21   1.4710   N21 180.0000    :1  180.0000
              H21 C22   1.0990   C21 109.8381   N21 280.0000
              H22 C22   1.0990   C21 109.8381   N21  20.4038
              H23 C22   1.0990   C21 109.8381   N21 140.4038
              0 1
          '''

        atom_name = [
          'N3', 'C2', 'C1', 'H11', 'H12', 'H13'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def monomer_b_carbon_z_matrix(self):

        zmatrix = '''\
                C21  :1   DISTANCE  :2  ANGLE    :3  DIHEDRAL
                X21 N21   1.0000    :1  90.0000    :2    0.0000
                N21 C21   1.1854   X21  90.0000    :1  180.0000
                X22 C21   1.0000   N21  90.0000   X21   0.00000
                C22 C21   1.4710   N21 180.0000    :1  180.0000
                H21 C22   1.0990   C21 109.8381   N21 280.0000
                H22 C22   1.0990   C21 109.8381   N21  20.4038
                H23 C22   1.0990   C21 109.8381   N21 140.4038
                0 1
            '''

        atom_name = [
          'C2', 'N3', 'C1', 'H11', 'H12', 'H13'
        ]

        return textwrap.dedent(zmatrix), atom_name

