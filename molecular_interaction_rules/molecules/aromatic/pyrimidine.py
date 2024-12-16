#!/usr/bin/env python3
#
# Lennard-Jones-Drill-2: Pyrimidine
# ---------------------------------

# Imports
# -------
import textwrap

class Pyrimidine(object):

    def __init__(self):

        self.resi_name = 'PYRM'

    def get_monomer_a_species(self):

        '''

        Get the Monomer A Species

        '''

        monomer_a_species = {
            'RC1': self.monomer_a_aromatic_zmatrix(),
            'N1': self.monomer_a_nitrogen_zmatrix()
        }

        return monomer_a_species

    def get_monomer_b_species(self):

        monomer_b_species = {
            'RC1': self.monomer_b_aromatic_zmatrix(),
        }

        return monomer_b_species

    def monomer_a_aromatic_zmatrix(self):

        zmatrix = '''\
            X11
            C11  X11  1.2940
            C12  C11  1.4026   X11   60.0000
            N11  C12  1.3508   C11  122.1105   X11    0.0000
            C13  N11  1.3497   C12  115.6860   C11    0.0000
            N12  C13  1.3497   N11  127.4170   C12    0.0000
            C14  C11  1.4026   C12  116.9770   N11   -0.0000
            H11  C13  1.0932   N11  116.2910   C12  180.0000
            H12  C14  1.0946   C11  116.9770   C12  180.0000
            H13  C12  1.0946   N11  116.4660   C13  180.0000
            H14  C11  1.0926   C12  121.5110   N11  180.0000
            0 1
        '''

        atom_name = [
            'C5', 'C6', 'N1', 'C2', 'N3', 'C4', 'H2', 'H4', 'H6', 'H5'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def monomer_a_nitrogen_zmatrix(self):

        zmatrix = '''\
            N11
            C12  N11  1.3508
            C11  C12  1.4026   N11  122.1105
            C13  N11  1.3497   C12  115.6860   C11    0.0000
            N12  C13  1.3497   N11  127.4170   C12    0.0000
            C14  C11  1.4026   C12  116.9770   N11   -0.0000
            H11  C13  1.0932   N11  116.2910   C12  180.0000
            H12  C14  1.0946   C11  116.9770   C12  180.0000
            H13  C12  1.0946   N11  116.4660   C13  180.0000
            H14  C11  1.0926   C12  121.5110   N11  180.0000
            X11  N11  1.0000   C12   90.0000   C11  180.0000
            0 1
        '''

        atom_name = [
          'N1', 'C6', 'C5', 'C2', 'N3', 'C4', 'H2', 'H4', 'H6', 'H5'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def monomer_b_aromatic_zmatrix(self):


      zmatrix = '''\
           X21   :1  DISTANCE   :2  ANGLE     :3   90.0000
            C21  X21  1.2940     :1  90.0000    :2   0.0000
            C22  C21  1.4026   X21   60.0000    :1   DIHEDRAL
            N21  C22  1.3508   C21  122.1105   X21    0.0000
            C23  N21  1.3497   C22  115.6860   C21    0.0000
            N22  C23  1.3497   N21  127.4170   C22    0.0000
            C24  C21  1.4026   C22  116.9770   N21   -0.0000
            H21  C23  1.0932   N21  116.2910   C22  180.0000
            H22  C24  1.0946   C21  116.9770   C22  180.0000
            H23  C22  1.0946   N21  116.4660   C23  180.0000
            H24  C21  1.0926   C22  121.5110   N21  180.0000
            0 1
        '''

      atom_name = [
        'C5', 'C6', 'N1', 'C2', 'N3', 'C4', 'H2', 'H4', 'H6', 'H5'
      ]

      return textwrap.dedent(zmatrix), atom_name

