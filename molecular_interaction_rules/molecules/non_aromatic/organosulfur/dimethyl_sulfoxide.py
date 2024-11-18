#!/usr/bin/env python3
#
# Lennard-Jones-Drill-2: DimethylSulfoxide
# ----------------------------------------

# Imports
# -------

import textwrap

class DimethylSulfoxide(object):

    def __init__(self):

        self.resi_name = 'DMSO'

    def get_monomer_a_species(self):

        '''

        Get the Monomer A Species

        '''

        monomer_a_species = {
            'S1': self.get_sulphur_zmatrix(),
            'O1': self.get_monomer_a_oxygen_zmatrix()
        }

        return monomer_a_species

    def get_monomer_a_oxygen_zmatrix(self):

        zmatrix = '''\
            O11
            S11  O11  1.5447
            C11  S11  1.8243  O11  130.0000
            C12  S11  1.8243  C11   95.8631  O11  180.0000
            H11  C12  1.1012  S11  109.5528  C11  -64.0096
            H12  C12  1.1005  S11  108.0603  C11   58.0274
            H13  C12  1.0989  S11  106.6409  C11  176.5567
            H14  C11  1.1012  S11  109.5528  C12   64.0096
            H15  C11  1.0989  S11  106.6409  C12 -176.5567
            H16  C11  1.1005  S11  108.0603  C12  -58.0274
            X11  S11  1.0000  C11  135.0000  O11  90.0000
            0 1
        '''

        atom_name = [
          'O1', 'S2', 'C3', 'C7', 'H8', 'H9', 'H10', 'H4', 'H5', 'H6'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_sulphur_zmatrix(self):

        zmatrix = '''\
            S11
            C11  S11  1.8243
            C12  S11  1.8243  C11   95.8631
            H11  C12  1.1012  S11  109.5528  C11  -64.0096
            H12  C12  1.1005  S11  108.0603  C11   58.0274
            H13  C12  1.0989  S11  106.6409  C11  176.5567
            O11  S11  1.5447  C12  105.9417  H11 -172.4598
            H14  C11  1.1012  S11  109.5528  C12   64.0096
            H15  C11  1.0989  S11  106.6409  C12 -176.5567
            H16  C11  1.1005  S11  108.0603  C12  -58.0274
            X11  S11  1.0000  C11  135.0000  O11  90.0000
            0 1
        '''

        atom_name = [
            'S2', 'C3', 'C7', 'H8', 'H9', 'H10', 'O1', 'H4', 'H5', 'H6'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_b_species(self):

        '''

        Get the Monomer B Species

        '''

        monomer_b_species = {
            'O1': self.get_monomer_b_oxygen_zmatrix()
        }

        return monomer_b_species

    def get_monomer_b_oxygen_zmatrix(self):

        zmatrix = '''\
            O21   :1  DISTANCE   :2 ANGLE    :3    DIHEDRAL
            S21  O21  1.5447     :1  180.0000    :2   180.0000
            C21  S21  1.8243  O21  130.0000      :1     0.0000
            C22  S21  1.8243  C21   95.8631  O21  180.0000
            H21  C22  1.1012  S21  109.5528  C21  -64.0096
            H22  C22  1.1005  S21  108.0603  C21   58.0274
            H23  C22  1.0989  S21  106.6409  C21  176.5567
            H24  C21  1.1012  S21  109.5528  C22   64.0096
            H25  C21  1.0989  S21  106.6409  C22 -176.5567
            H26  C21  1.1005  S21  108.0603  C22  -58.0274
            X21  S21  1.0000  C21  135.0000  O21  90.0000
            0 1
        '''

        atom_name = [
          'O1', 'S2', 'C3', 'C7', 'H8', 'H9', 'H10', 'H4', 'H5', 'H6'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_b_sulphur_zmatrix(self):

        zmatrix = '''\
            S21
            C21  S21  1.8243
            C22  S21  1.8243  C21   95.8631
            H21  C22  1.1012  S21  109.5528  C21  -64.0096
            H22  C22  1.1005  S21  108.0603  C21   58.0274
            H23  C22  1.0989  S21  106.6409  C21  176.5567
            O21  S21  1.5447  C22  105.9417  H21 -172.4598
            H24  C21  1.1012  S21  109.5528  C22   64.0096
            H25  C21  1.0989  S21  106.6409  C22 -176.5567
            H26  C21  1.1005  S21  108.0603  C22  -58.0274
            X21  S21  1.0000  C21  135.0000  O21  90.0000
            0 1
        '''

        atom_name = [
          'S2', 'C3', 'C7', 'H8', 'H9', 'H10', 'O1', 'H4', 'H5', 'H6'
        ]

        return textwrap.dedent(zmatrix), atom_name
