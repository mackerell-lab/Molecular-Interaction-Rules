#!/usr/bin/env python3
#
# Lennard-Jones-Drill-2: Hydrazine
# --------------------------------

# Imports
# -------
import textwrap

class Hydrazine(object):

    def __init__(self):

        self.resi_name = 'HDZN'

    def get_monomer_a_species(self):

        '''

        Get the Monomer A Species

        '''

        monomer_a_species = {
            'H1': self.get_monomer_a_nitrogen_zmatrix()
        }

        return monomer_a_species

    def get_monomer_b_species(self):

        monomer_b_species = {
            'H1': self.get_monomer_b_nitrogen_zmatrix()
        }

        return monomer_b_species

    def get_monomer_a_nitrogen_zmatrix(self):

        zmatrix = '''\
            H11
            N12 H11   1.0187
            N11 N12   1.4446   H11 111.4203
            H12  N12  1.0187  N11  111.4203  H11  116.4207
            H13  N11  1.0187  N12  111.4203  H11  -89.9827
            H14  N11  1.0187  N12  111.4203  H11  153.5965
            X11  N11  1.0000  N12   90.0000  H11  180.0000
            0 1
        '''

        atom_name = [
          'H11', 'N1', 'N2',  'H12', 'H21', 'H22'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_b_nitrogen_zmatrix(self):

        zmatrix = '''\
            H21   :1   DISTANCE     :2  ANGLE   :3   DIHEDRAL
            X21  H21    1.0000      :1  90.0000   :2    0.0000
            N22  H21   1.0219   X21   90.0000     :2  180.0000
            N21  N22   1.4446   H21 111.4203      :1  180.0000
            H22  N22  1.0187    N21  111.4203  H21  116.4207
            H23  N21  1.0187    N22  111.4203  H21  -89.9827
            H24  N21  1.0187    N22  111.4203  H21  153.5965
            X21  N21  1.0000    N22   90.0000  H21  180.0000
            0 1
        '''

        atom_name = [
          'H11', 'N1', 'N2', 'H12', 'H21', 'H22'
        ]

        return textwrap.dedent(zmatrix), atom_name


