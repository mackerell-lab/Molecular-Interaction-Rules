#!/usr/bin/env python3
#
# Lennard-Jones-Drill-2: Azetidinone
# ----------------------------------

# Imports
# -------

import textwrap

class Azetidinone(object):

    def __init__(self):

        self.resi_name = 'AZDO'

    def get_monomer_a_species(self):

        monomer_a_species = {
            'H1': self.get_monomer_a_nitrogen_zmatrix(),
            'O1': self.get_monomer_a_oxygen_zmatrix(),
        }

        return monomer_a_species

    def get_monomer_b_species(self):

      '''

      Get the Monomer B Species

      '''

      monomer_b_species = {
        'O1': self.get_monomer_b_oxygen_zmatrix(),
        'H1': self.get_monomer_b_hydrogen_zmatrix()
      }

      return monomer_b_species

    def get_monomer_a_nitrogen_zmatrix(self):

        zmatrix = '''\
          H11
          N11  H11  1.0158
          C12  N11  1.3624  H11  130.7810
          C11  C12  1.5314  N11   91.3607  H11 -161.6569
          C13  C11  1.5484  C12   85.4333  N11   -5.5311
          O11  C12  1.1855  C11  135.8730  C13  173.2256
          H12  C13  1.0837  C11  114.3650  C12 -110.0105
          H13  C13  1.0821  C11  116.1886  C12  120.6051
          H14  C11  1.0815  C13  116.2706  N11  120.1111
          H15  C11  1.0819  C13  115.3314  N11 -109.1013
          X11  H11  1.0000  N11   90.0000  C12  180.0000
          0 1
        '''

        atom_name = [
          'H1', 'N1', 'C2', 'C3', 'C4', 'O2', 'H41', 'H42', 'H32', 'H31'
        ]

        return zmatrix, atom_name

    def get_monomer_a_oxygen_zmatrix(self):

        zmatrix = '''\
          O11
          C12  O11  1.1855
          C11  C12  1.5314  O11  135.8730
          C13  C11  1.5484  C12   85.4333  O11  173.2256
          N11  C12  1.3624  C11   91.3607  C13   -5.5311
          H11  N11  1.0158  C12  130.7810  C11 -161.6569
          H12  C13  1.0837  C11  114.3650  C12 -110.0105
          H13  C13  1.0821  C11  116.1886  C12  120.6051
          H14  C11  1.0815  C13  116.2706  N11  120.1111
          H15  C11  1.0819  C13  115.3314  N11 -109.1013
          X11  O11  1.0000  C12   90.0000  C11  180.0000
          0 1
        '''

        atom_name = [
          'O2', 'C2', 'C3', 'C4', 'N1', 'H1', 'H41', 'H42', 'H32', 'H31'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_b_hydrogen_zmatrix(self):

        zmatrix = '''\
            H21   :1    DISTANCE   :2 ANGLE     :3 DIHEDRAL
            X21  H21  1.0000   :1  90.0000     :2    0.0000
            N21  H21  1.0158  X21  90.0000     :2  180.0000
            C22  N21  1.3624  H21  130.7810    :1  180.0000
            C21  C22  1.5314  N21   91.3607  H21 -161.6569
            C23  C21  1.5484  C22   85.4333  N21   -5.5311
            O21  C22  1.1855  C21  135.8730  C23  173.2256
            H22  C23  1.0837  C21  114.3650  C22 -110.0105
            H23  C23  1.0821  C21  116.1886  C22  120.6051
            H24  C21  1.0815  C23  116.2706  N21  120.1111
            H25  C21  1.0819  C23  115.3314  N21 -109.1013
            0 1
        '''

        atom_name = [
          'H1', 'N1', 'C2', 'C3', 'C4', 'O2', 'H41', 'H42', 'H32', 'H31'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_b_oxygen_zmatrix(self):

        zmatrix = '''\
            O21   :1  DISTANCE  :2  ANGLE     :3   90.0000
            C22  O21  1.1855    :1  180.0000  :2   DIHEDRAL
            C21  C22  1.5314  O21  135.8730   :1    0.0000
            C23  C21  1.5484  C22   85.4333  O21  173.2256
            N21  C22  1.3624  C21   91.3607  C23   -5.5311
            H21  N21  1.0158  C22  130.7810  C21 -161.6569
            H22  C23  1.0837  C21  114.3650  C22 -110.0105
            H23  C23  1.0821  C21  116.1886  C22  120.6051
            H24  C21  1.0815  C23  116.2706  N21  120.1111
            H25  C21  1.0819  C23  115.3314  N21 -109.1013
            X21  O21  1.0000  C22   90.0000  C21  180.0000
            0 1
          '''

        atom_name = [
          'O2', 'C2', 'C3', 'C4', 'N1', 'H1', 'H41', 'H42', 'H32', 'H31'
        ]

        return textwrap.dedent(zmatrix), atom_name

