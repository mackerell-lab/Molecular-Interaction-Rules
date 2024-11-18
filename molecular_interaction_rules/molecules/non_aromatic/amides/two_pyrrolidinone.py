#!/usr/bin/env python3
#
# Lennard-Jones-Drill-2: Two Pyrrolidinone
# ----------------------------------------

# Imports
# -------

import textwrap

class TwoPyrrolidinone(object):

    def __init__(self):

        self.resi_name = '2PD0'

    def get_monomer_a_species(self):

        '''

        Get the Monomer A Species

        '''

        monomer_a_species = {
            'H1': self.get_monomer_a_nitrogen_zmatrix()
        }

        return monomer_a_species

    def get_monomer_a_nitrogen_zmatrix(self):

        zmatrix = '''\
            H16
            N11   H16   1.0146
            C11   N11   1.3801   H16  119.0662
            C12   C11   1.5203   N11  107.0109    H16 -169.7098
            C13   C12   1.5334   C11  104.2838    N11  -12.6895
            H11   C13   1.0824   C12  113.4706    C11  147.0441
            H12   C13   1.0842   C12  109.9870    C11  -91.9544
            C14   N11   1.4472   C11  114.7839    C12   -6.0332
            H12   C14   1.0876   N11  111.5487    C11  -97.3756
            H13   C14   1.0827   N11  111.1615    C11  142.0328
            H14   C12   1.0815   C11  110.6249    N11 -136.2340
            H15   C12   1.0862   C11  107.5780    N11  106.4460
            O11   C11   1.2303   C12  125.8982    C13  168.1336
            X11   H16   1.0000   N11   90.0000    C11    0.0000
            0 1
        '''

        atom_name = [
            'H1', 'N1', 'C2', 'C3', 'C4', 'H41', 'H42', 'C5', 'H51', 'H52', 'C31', 'C32', 'O2'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_b_species(self):

        '''

        Get the Monomer B Species

        '''

        monomer_b_species = {
            'H1': self.get_monomer_b_nitrogen_zmatrix()
        }

        return monomer_b_species

    def get_monomer_b_nitrogen_zmatrix(self):

        zmatrix = '''\
              H26    :1   DISTANCE  :2   ANGLE   :3   DIHEDRAL
              X21   H26    1.0000    :1   90.0000   :2    0.0000
              N21   H26   1.0146   X21  90.0000    :1  180.0000
              C21   N21   1.3801   H26  119.5782   :1  180.0000
              C22   C21   1.5203   N21  107.0109    H26 -169.7098
              C23   C22   1.5334   C21  104.2838    N21  -12.6895
              H21   C23   1.0824   C22  113.4706    C21  147.0441
              H22   C23   1.0842   C22  109.9870    C21  -91.9544
              C24   N21   1.4472   C21  114.7839    C22   -6.0332
              H22   C24   1.0876   N21  111.5487    C21  -97.3756
              H23   C24   1.0827   N21  111.1615    C21  142.0328
              H24   C22   1.0815   C21  110.6249    N21 -136.2340
              H25   C22   1.0862   C21  107.5780    N21  106.4460
              O21   C21   1.1958   C22  126.9545    C23  168.1336
              0 1
          '''

        atom_name = [
          'H1', 'N1', 'C2', 'C3', 'C4', 'H41', 'H42', 'C5', 'H51', 'H52', 'C31', 'C32', 'O2'
        ]

        return textwrap.dedent(zmatrix), atom_name



