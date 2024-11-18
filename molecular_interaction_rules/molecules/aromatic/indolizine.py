#!/usr/bin/env python3
#
# Lennard-Jones-Drill-2: Indolizine
# ---------------------------------

# Imports
# -------
import textwrap

class Indolizine(object):

    def __init__(self):

        self.resi_name = 'indz'

    def get_monomer_a_species(self):

        monomer_a_species = {
            'RC1': self.monomer_a_aromatic_zmatrix(),
            # 'N1': self.monomer_a_nitrogen_zmatrix()
        }

        return monomer_a_species

    def get_monomer_b_species(self):

        monomer_b_species = {
            'RC1': self.monomer_b_aromatic_zmatrix()
        }

        return monomer_b_species

    def monomer_a_aromatic_zmatrix(self):

        zmatrix = '''\
          X11
          N11  X11 1.1940
          C11  N11  1.3828  X11   60.0000
          C12  C11  1.3803  N11  119.7565   X11    0.0000
          C13  C12  1.4302  C11  120.6364   N11    0.0000
          C14  C13  1.3904  C12  119.4765   C11   -0.0000
          C15  N11  1.4230  C11  121.5132   C12   -0.0000
          C16  C15  1.4133  N11  106.7973   C11 -180.0000
          C17  C16  1.4148  C15  107.3590   C14 -180.0000
          C18  N11  1.3770  C11  128.9684   C12  180.0000
          H11  C13  1.0927  C12  120.1455   C11  180.0000
          H12  C14  1.0934  C13  121.1538   C12 -180.0000
          H13  C16  1.0889  C15  125.1425   N11 -180.0000
          H14  C17  1.0899  C16  126.6314   C15 -180.0000
          H15  C18  1.0884  N11  120.8740   C11   -0.0000
          H16  C12  1.0922  C11  118.3940   N11 -180.0000
          H17  C11  1.0918  C12  124.0410   C13 -180.0000
          X11  N11  1.0000  C11   90.0000   C12  180.0000
          0 1
        '''

        atom_name = [
          'N4', 'C5', 'C6', 'C7', 'C8', 'C9', 'C1', 'C2', 'C3', 'H7', 'H8', 'H4', 'H2', 'H3', 'H6', 'H5',
        ]

        return textwrap.dedent(zmatrix), atom_name

    def monomer_b_aromatic_zmatrix(self):

        zmatrix = '''\
            X21   :1  DISTANCE  :2   ANGLE     :3   90.0000
            N21  X21 1.1940     :1   90.0000   :2    0.0000
            C21  N21  1.3828  X21   60.0000    :1   DIHEDRAL
            C22  C21  1.3803  N21  119.7565   X21    0.0000
            C23  C22  1.4302  C21  120.6364   N21    0.0000
            C24  C23  1.3904  C22  119.4765   C21   -0.0000
            C25  N21  1.4230  C21  121.5132   C22   -0.0000
            C26  C25  1.4133  N21  106.7973   C21 -180.0000
            C27  C26  1.4148  C25  107.3590   C24 -180.0000
            C28  N21  1.3770  C21  128.9684   C22  180.0000
            H21  C23  1.0927  C22  120.1455   C21  180.0000
            H22  C24  1.0934  C23  121.1538   C22 -180.0000
            H23  C26  1.0889  C25  125.1425   N21 -180.0000
            H24  C27  1.0899  C26  126.6314   C25 -180.0000
            H25  C28  1.0884  N21  120.8740   C21   -0.0000
            H26  C22  1.0922  C21  118.3940   N21 -180.0000
            H27  C21  1.0918  C22  124.0410   C23 -180.0000
            0 1
          '''

        atom_name = [
          'N4', 'C5', 'C6', 'C7', 'C8', 'C9', 'C1', 'C2', 'C3', 'H7', 'H8', 'H4', 'H2', 'H3', 'H6', 'H5',
        ]

        return textwrap.dedent(zmatrix), atom_name

    def monomer_a_nitrogen_zmatrix(self):

        zmatrix = '''\
            N11
            C11  N11  1.3828
            C12  C11  1.3803  N11  119.7565
            C13  C12  1.4302  C11  120.6364   N11    0.0000
            C14  C13  1.3904  C12  119.4765   C11   -0.0000
            C15  N11  1.4230  C11  121.5132   C12   -0.0000
            C16  C15  1.4133  N11  106.7973   C11 -180.0000
            C17  C16  1.4148  C15  107.3590   C14 -180.0000
            C18  N11  1.3770  C11  128.9684   C12  180.0000
            H11  C13  1.0927  C12  120.1455   C11  180.0000
            H12  C14  1.0934  C13  121.1538   C12 -180.0000
            H13  C16  1.0889  C15  125.1425   N11 -180.0000
            H14  C17  1.0899  C16  126.6314   C15 -180.0000
            H15  C18  1.0884  N11  120.8740   C11   -0.0000
            H16  C12  1.0922  C11  118.3940   N11 -180.0000
            H17  C11  1.0918  C12  124.0410   C13 -180.0000
            X11  N11  1.0000  C11   90.0000   C12  180.0000
            0 1
        '''

        atom_name = [
            'N4', 'C5', 'C6', 'C7', 'C8', 'C9', 'C1', 'C2', 'C3', 'H7', 'H8', 'H4', 'H2', 'H3', 'H6', 'H5',
        ]

        return textwrap.dedent(zmatrix), atom_name

    def monomer_b_pi_stack_zmatrix(self):


        zmatrix = '''\
        '''

        atom_name = [
        ],

        return textwrap.dedent(zmatrix), atom_name

