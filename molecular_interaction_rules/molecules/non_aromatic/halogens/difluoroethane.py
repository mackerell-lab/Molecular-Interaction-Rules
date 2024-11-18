#!/usr/bin/env python3
#
# Lennard-Jones-Drill-2: Difluoroethane
# -------------------------------------

# Imports
# -------

import textwrap

class Difluoroethane(object):

    def __init__(self):

        self.resi_name = 'DFET'

    def get_monomer_a_species(self):

        monomer_a_species = {
            'F1': self.get_monomer_a_fluoro_zmatrix(),
            'H1': self.get_monomer_a_hydrogen_zmatrix()
        }

        return monomer_a_species

    def get_monomer_a_fluoro_zmatrix(self):

        zmatrix = '''\
            F11
            C11   F11   1.3863
            C12   C11   1.5054   F11  110.1466
            H11   C12   1.0996   C11  109.4854    F11  -58.6258
            H12   C12   1.0996   C11  109.4854    F11 -178.6928
            H13   C12   1.0996   C11  109.4854    F11   61.4412
            H14   C11   1.0996   C12  115.1421    H11 -180.0000
            F12   C11   1.3863   C12  110.1466    H11   58.6258
            X11   F11   1.0000   C11   90.0000    C12  180.0000
            0 1
        '''

        atom_name = [
            'F22', 'C2', 'C1', 'H11', 'H12', 'H13', 'H21', 'F21'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_a_hydrogen_zmatrix(self):

        zmatrix = '''\
              H11
              C11  H11  1.0996
              C12  C11  1.5054  H11  109.4854
              F11  C11  1.3863  C12  110.1466   H11  115.4412
              H12  C12  1.0996  C11  109.4854   F11  -58.6258
              H13  C12  1.0996  C11  109.4854   F11 -178.6928
              H14  C12  1.0996  C11  109.4854   F11   61.4412
              F12  C11  1.3863  C12  110.1466   H12   58.6258
              X11  H11  1.0000  C11   90.0000   C12  180.0000
              0 1
          '''

        atom_name = [
          'H21', 'C2', 'C1', 'F22', 'H11', 'H12', 'H13', 'H21', 'F21'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_b_species(self):

        monomer_b_species = {
          'F1': self.get_monomer_a_fluoro_zmatrix(),
          'H1': self.get_monomer_b_hydrogen_zmatrix()
        }

        return monomer_b_species

    def get_monomer_b_fluoro_zmatrix(self):

        zmatrix = '''\
            F21   :1 DISTANCE    :2  ANGLE      :3  DIHEDRAL
            X21   F21 1.0000      :1  90.0000     :2   0.0000
            C21   F21   1.3863   X21  90.0000    :1  180.0000
            C22   C21   1.5054   F21  110.1466   :1  180.0000
            H21   C22   1.0996   C21  109.4854    F21  -58.6258
            H22   C22   1.0996   C21  109.4854    F21 -178.6928
            H23   C22   1.0996   C21  109.4854    F21   61.4412
            H24   C21   1.0996   C22  115.1421    H21 -180.0000
            F22   C21   1.3863   C22  110.1466    H21   58.6258
            0 1
        '''

        atom_name = [
          'F22', 'C2', 'C1', 'H11', 'H12', 'H13', 'H21', 'F21'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_b_hydrogen_zmatrix(self):

        zmatrix = '''\
                H21   :1  DISTANCE    :2  ANGLE   :3  DIHEDRAL
                X21  H21  1.0000      :1  90.0000   :2   0.0000
                C21  H21  1.0996  X21  90.0000      :1  180.0000
                C22  C21  1.5054  H21  109.4854     :1  180.0000
                F21  C21  1.3863  C22  110.1466   H21  115.4412
                H22  C22  1.0996  C21  109.4854   F21  -58.6258
                H23  C22  1.0996  C21  109.4854   F21 -178.6928
                H24  C22  1.0996  C21  109.4854   F21   61.4412
                F22  C21  1.3863  C22  110.1466   H22   58.6258
                0 1
            '''

        atom_name = [
             'H21', 'C2', 'C1', 'F22', 'H11', 'H12', 'H13', 'H21', 'F21'
        ]

        return textwrap.dedent(zmatrix), atom_name

