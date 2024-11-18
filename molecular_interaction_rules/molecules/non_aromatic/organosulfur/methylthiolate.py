#!/usr/bin/env python3
#
# Lennard-Jones-Drill-2: MethylThiolate
# -------------------------------------

# Imports
# -------

import textwrap

class MethylThiolate(object):

    def __init__(self):

        self.resi_name = 'MES1'

    def get_monomer_a_species(self):

        '''

        Get the Monomer A Species

        '''

        monomer_a_species = {
          'H1': self.get_monomer_a_hydrogen_zmatrix(),
          'S1': self.get_monomer_a_sulphur_zmatrix()
        }

        return monomer_a_species

    def get_monomer_a_hydrogen_zmatrix(self):

        zmatrix = '''\
          H11
          C11  H11 1.0990
          S11  C11 1.8465 H11  106.1319
          H12  C11 1.1067 S11  111.1621  H11  -118.2020
          H13  C11 1.1067 S11  111.1621  H11   118.2020
          X11  H11 1.0000 C11   90.0000  S11   180.0000
          -1 1
        '''

        atom_name = [
          'H11', 'C1', 'S', 'H12', 'H13'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_a_sulphur_zmatrix(self):

        zmatrix = '''\
          S11
          C11  S11 1.8465
          H11  C11 1.1067 S11  106.1319
          H12  C11 1.1067 S11  111.1621  H11 -118.2020
          H13  C11 1.1067 S11  111.1621  H11  118.2020
          X11  S11 1.0000 C11   90.0000  H11  180.0000
          -1 1
        '''

        atom_name = [
          'S', 'C1', 'H11', 'H12', 'H13'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_b_species(self):

        '''

        Get the Monomer B Species

        '''

        monomer_b_species = {
            'H1': self.get_monomer_b_hydrogen_zmatrix(),
            'S1': self.get_monomer_b_sulphur_zmatrix()
        }

        return monomer_b_species

    def get_monomer_b_hydrogen_zmatrix(self):

        zmatrix = '''\
            H21   :1    DISTANCE     :2 ANGLE    :3 DIHEDRAL
            X21   H21   1.0000     :1  90.0000    :2    0.0000
            C21  H21 1.0990   X21  90.0000          :2  180.0000
            S21  C21 1.8465 H21  106.1319           :1  180.0000
            H22  C21 1.1067 S21  111.1621  H21  -118.2020
            H23  C21 1.1067 S21  111.1621  H21   118.2020
            -1 1
          '''

        atom_name = [
          'H11', 'C1', 'S', 'H12', 'H13'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_b_sulphur_zmatrix(self):

        zmatrix = '''\
          S21   :1  DISTANCE   :2  ANGLE    :3    DIHEDRAL
          C21  S21 1.8465      :1  180.0000    :2   180.0000
          H21  C21 1.1067 S21  106.1319        :1     0.0000
          H22  C21 1.1067 S21  111.1621  H21 -118.2020
          H23  C21 1.1067 S21  111.1621  H21  118.2020
          X21  S21 1.0000 C21   90.0000  H21  180.0000
          -1 1
        '''

        atom_name = [
          'S', 'C1', 'H11', 'H12', 'H13'
        ]

        return textwrap.dedent(zmatrix), atom_name

