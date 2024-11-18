#!/usr/bin/env python3
#
# Lennard-Jones-Drill-2: Tetrahydropyran
# --------------------------------------

# Imports
# -------

import textwrap

class Tetrahydropyran(object):

    def __init__(self):

        self.resi_name = 'THP'

    def get_monomer_a_species(self):

        '''

        Get the Monomer A Species

        '''

        monomer_a_species = {
            'O1': self.get_monomer_a_oxygen_zmatrix()
        }

        return monomer_a_species

    def get_monomer_a_oxygen_zmatrix(self):

        zmatrix = '''\
          O11
          C11  O11  1.4342
          C12  C11  1.5423  O11  110.2350
          C13  C12  1.5423  C11  109.3871  O11   31.6498
          H11  C13  1.0994  C12  109.7015  C11  -85.8318
          H12  C13  1.0994  C12  110.7313  C11  155.8732
          C14  C13  1.5345  C12  108.7313  C11   33.3583
          H13  C14  1.0994  C13  111.8031  C12  174.9336
          H14  C14  1.0994  C13  109.2669  C12   56.7576
          C15  C14  1.5432  C13  108.7987  C12  -64.1404
          H15  C15  1.0994  O11  104.9354  C11  159.5652
          H16  C15  1.0994  O11  109.2669  C11  -84.1985
          H17  C12  1.0994  C11  110.7313  O11  -89.6040
          H18  C12  1.0994  C11  109.2669  O11  152.5557
          H19  C11  1.0994  C12  111.8902  C13  149.6187
          H20  C11  1.1078  C12  111.8902  C13  -89.8718
          X11  O11  1.0000  C11   90.0000  C12  180.0000
          0 1
        '''

        atom_name = [
          'O1', 'C2', 'C3', 'H31', 'H32', 'C4', 'H41', 'H42', 'C5', 'H51', 'H52', 'H21', 'H22', 'H11', 'H12'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_b_species(self):

        monomer_b_species = {
            'O1': self.get_monomer_b_oxygen_zmatrix()
        }

        return monomer_b_species

    def get_monomer_b_oxygen_zmatrix(self):

        zmatrix = '''\
            O21   :1  DISTANCE  :2   ANGLE       :3  DIHEDRAL
            C21  O21  1.4342    :1   124.4855    :2   0.00000
            C22  O21  1.5423    :1   124.4855    :2  180.00000
            C23  C22  1.5423  C21  109.3871  O21   31.6498
            H21  C23  1.0994  C22  109.7015  C21  -85.8318
            H22  C23  1.0994  C22  110.7313  C21  155.8732
            C24  C23  1.5345  C22  108.7313  C21   33.3583
            H23  C24  1.0994  C23  111.8031  C22  174.9336
            H24  C24  1.0994  C23  109.2669  C22   56.7576
            C25  C24  1.5432  C23  108.7987  C22  -64.1404
            H25  C25  1.0994  O21  104.9354  C21  159.5652
            H26  C25  1.0994  O21  109.2669  C21  -84.1985
            H27  C22  1.0994  C21  110.7313  O21  -89.6040
            H28  C22  1.0994  C21  109.2669  O21  152.5557
            H29  C21  1.0994  C22  111.8902  C23  149.6187
            H30  C21  1.1078  C22  111.8902  C23  -89.8718
            0 1
          '''

        atom_name = [
          'O1', 'C2', 'C3', 'H31', 'H32', 'C4', 'H41', 'H42', 'C5', 'H51', 'H52', 'H21', 'H22', 'H11', 'H12'
        ]

        return textwrap.dedent(zmatrix), atom_name
