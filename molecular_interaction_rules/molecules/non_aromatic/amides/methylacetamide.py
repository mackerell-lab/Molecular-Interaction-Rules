#!/usr/bin/env python3
#
# Lennard-Jones-Drill-2: MethylAcetamide
# --------------------------------------

# Imports
# -------

import textwrap

class MethylAcetamide(object):

    def __init__(self):

        self.resi_name = 'NMA'

    def get_monomer_a_species(self):

        '''

        Get the Monomer A Species

        '''

        monomer_a_species = {
            'N1': self.get_monomer_a_nitrogen_zmatrix()
        }

        return monomer_a_species

    def get_monomer_a_nitrogen_zmatrix(self):

        zmatrix = '''\
            H17
            N11  H17  1.0104
            C11  N11  1.3686   H17  116.5274
            C12  C11  1.5185   N11  116.1617   H17    0.0000
            H11  C12  1.0996   C11  108.5139   N11 -121.4890
            H12  C12  1.0996   C11  108.5134   N11  121.5323
            H13  C12  1.1001   C11  113.2240   N11    0.0216
            O11  C11  1.2378   C12  122.4208   H11   58.5138
            C13  N11  1.4579   C11  120.3664   C12  179.9797
            H14  C13  1.1001   N11  110.5107   C11  -59.8292
            H15  C13  1.1001   N11  110.5001   C11   59.7654
            H16  C13  1.0978   N11  108.6683   C11  179.9643
            X11  H17  1.0000   N11   90.0000   C11  180.0000
            0 1
        '''

        atom_name = [
          'HR3', 'N', 'C', 'CL', 'HL1', 'HL2', 'HL3', 'O', 'CR', 'HR1', 'HR2',
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_b_species(self):

        '''

        Get the Monomer B Species

        '''

        monomer_b_species = {
            'N1': self.get_monomer_b_nitrogen_zmatrix()
        }

        return monomer_b_species

    def get_monomer_b_nitrogen_zmatrix(self):

        zmatrix = '''\
            H27   :1  DISTANCE  :2  ANGLE     :3  DIHEDRAL
            X21  H27  1.0000   :1   90.0000    :2 0.0000
            N21  H27  1.0104   X21   90.0000   :1 180.0000
            C21  N21  1.3686   H27  116.5274   X21   0.0000
            C22  C21  1.5185   N21  116.1617   H27    0.0000
            H21  C22  1.0996   C21  108.5139   N21 -121.4890
            H22  C22  1.0996   C21  108.5134   N21  121.5323
            H23  C22  1.1001   C21  113.2240   N21    0.0216
            O21  C21  1.2378   C22  122.4208   H21   58.5138
            C23  N21  1.4579   C21  120.3664   C22  179.9797
            H24  C23  1.1001   N21  110.5107   C21  -59.8292
            H25  C23  1.1001   N21  110.5001   C21   59.7654
            H26  C23  1.0978   N21  108.6683   C21  179.9643
            X21  H27  1.0000   N21   90.0000   C21  180.0000
            0 1
          '''

        atom_name = [
          'HR3', 'N', 'C', 'CL', 'HL1', 'HL2', 'HL3', 'O', 'CR', 'HR1', 'HR2',
        ]

        return textwrap.dedent(zmatrix), atom_name

