#!/usr/bin/env python3
#
# Lennard-Jones-Drill-2: Thiophene
# --------------------------------

# Imports
# -------
import textwrap

class Thiophene(object):

    def __init__(self):

        self.resi_name = 'thip'

    def get_monomer_a_species(self):

        '''

        Get the Monomer A Species

        '''

        monomer_a_species = {
            'RC1': self.monomer_a_aromatic_zmatrix(),
            'S1': self.monomer_a_sulfur_zmatrix()
        }

        return monomer_a_species

    def get_monomer_b_species(self):

        monomer_b_species = {
          'RC1': self.monomer_b_aromatic_zmatrix(),
          'S1': self.monomer_b_sulfur_zmatrix()
        }

        return monomer_b_species

    def monomer_a_aromatic_zmatrix(self):

        zmatrix = '''\
            X11
            C11  X11  1.1000
            S11  C11  1.7300  X11   60.0000
            C12  S11  1.7300  C11   92.0682  X11    0.0000
            C13  C12  1.3916  S11  111.5010  C11   -0.0000
            C14  C12  2.3430  S11   77.2544  C11    0.0000
            H11  C13  1.0920  C12  122.9446  S11 -180.0000
            H12  C14  1.0920  C11  122.9447  S11 -180.0000
            H13  C12  1.0900  S11  119.9886  C11  180.0000
            H14  C11  1.0900  S11  119.9886  C12 -180.0000
            0 1
        '''

        atom_name = [
          'C1', 'S5', 'C4', 'C3', 'C2', 'H3', 'H2', 'H4', 'H1'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def monomer_a_sulfur_zmatrix(self):

        zmatrix = '''\
              S11
              C11  S11  1.7300
              C12  S11  1.7300  C11   92.0682
              C13  C12  1.3916  S11  111.5010  C11   -0.0000
              C14  C12  2.3430  S11   77.2544  C11    0.0000
              H11  C13  1.0920  C12  122.9446  S11 -180.0000
              H12  C14  1.0920  C11  122.9447  S11 -180.0000
              H13  C12  1.0900  S11  119.9886  C11  180.0000
              H14  C11  1.0900  S11  119.9886  C12 -180.0000
              X11  S11  1.0000  C11   90.0000  C12    0.0000
              0 1
          '''

        atom_name = [
           'S5', 'C1', 'C4', 'C3', 'C2', 'H3', 'H2', 'H4', 'H1'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def monomer_b_aromatic_zmatrix(self):

        zmatrix = '''\
              X21   :1  DISTANCE  :2  ANGLE    :3   90.0000
              C21  X21  1.1000    :1  90.0000    :2  180.0000
              S21  C21  1.7300  X21   60.0000   :1   90.0000
              C22  S21  1.7300  C21   92.0682  X21  DIHEDRAL
              C23  C22  1.3916  S21  111.5010  C21   -0.0000
              C24  C22  2.3430  S21   77.2544  C21    0.0000
              H21  C23  1.0920  C22  122.9446  S21 -180.0000
              H22  C24  1.0920  C21  122.9447  S21 -180.0000
              H23  C22  1.0900  S21  119.9886  C21  180.0000
              H24  C21  1.0900  S21  119.9886  C22 -180.0000
              0 1
          '''

        atom_name = [
          'C1', 'S5', 'C4', 'C3', 'C2', 'H3', 'H2', 'H4', 'H1'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def monomer_b_sulfur_zmatrix(self):


        zmatrix = '''\
              S21   :1  DISTANCE    :2  ANGLE   :3   DIHEDRAL
              C21  S21  1.7300   :1  126.5480   :2     0.0000
              C22  S21  1.7300  C21   92.0682   :1     0.0000
              C23  C22  1.3916  S21  111.5010  C21   -0.0000
              C24  C22  2.3430  S21   77.2544  C21    0.0000
              H21  C23  1.0920  C22  122.9446  S21 -180.0000
              H22  C24  1.0920  C21  122.9447  S21 -180.0000
              H23  C22  1.0900  S21  119.9886  C21  180.0000
              H24  C21  1.0900  S21  119.9886  C22 -180.0000
              X21  S21  1.0000  C21   90.0000  C22    0.0000
              0 1
          '''

        atom_name = [
          'S5', 'C1', 'C4', 'C3', 'C2', 'H3', 'H2', 'H4', 'H1'
        ]

        return textwrap.dedent(zmatrix), atom_name



