#!/usr/bin/env python3
#
# Lennard-Jones-Drill-2: Fluoroethane
# -----------------------------------

# Imports
# -------

import textwrap

class Fluoroethane(object):

    def __init__(self):

        self.resi_name = 'FETH'

    def get_monomer_a_species(self):

        '''

        Get the Monomer A Species

        '''

        monomer_a_species = {
            'F1': self.get_monomer_a_fluoro_zmatrix(),
            'H1': self.get_monomer_a_carbon_zmatrix()
        }

        return monomer_a_species

    def get_monomer_a_fluoro_zmatrix(self):

        zmatrix = '''\
            F11
            C11 F11 1.4174
            C12 C11 1.5152 F11 109.5365
            H11 C12 1.0996 C11 110.5586 F11   60.3866
            H12 C12 1.0996 C11 110.5586 F11  -60.3866
            H13 C12 1.1020 C11 109.4385 F11 -180.0000
            H14 C11 1.1006 C12 112.0835 H11  -57.8794
            H15 C11 1.1006 C12 112.0835 H11  178.6526
            X11 F11 1.0000 C11  90.0000 C12  180.0000
            0 1
        '''

        atom_name = [
            'F21', 'C2', 'C1', 'H11', 'H12', 'H13', 'H22', 'H23'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_a_carbon_zmatrix(self):

        zmatrix = '''\
              H15
              C11 H15 1.1006
              F11 C11 1.4174 H15 112.0835
              C12 C11 1.5152 F11 109.5365 H15   100.0000
              H11 C12 1.0996 C11 110.5586 F11   60.3866
              H12 C12 1.0996 C11 110.5586 F11  -60.3866
              H13 C12 1.1020 C11 109.4385 F11 -180.0000
              H14 C11 1.1006 C12 112.0835 H11  -57.8794
              X11 H15 1.0000 C11  90.0000 F11    0.0000
              0 1
          '''

        atom_name = [
          'H22', 'C2', 'F21','C1', 'H11', 'H12', 'H13',  'H23'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_b_species(self):

        monomer_b_species = {
            'F1': self.get_monomer_b_fluoroethane_zmatrix(),
            'H1': self.get_monomer_b_carbon_zmatrix()
        }

        return monomer_b_species

    def get_monomer_b_fluoroethane_zmatrix(self):

        zmatrix = '''\
            F21  :1 DISTANCE  :2 ANGLE  :3  DIHEDRAL
            X21 F21 1.0000  :1  90.0000   :2   0.0000
            C21 F21 1.4174 X21  90.0000   :1  180.0000
            C22 C21 1.5152 F21 109.5365   :1  180.0000
            H21 C22 1.0996 C21 110.5586 F21   60.3866
            H22 C22 1.0996 C21 110.5586 F21  -60.3866
            H23 C22 1.1020 C21 109.4385 F21 -180.0000
            H24 C21 1.1006 C22 112.0835 H21  -57.8794
            H25 C21 1.1006 C22 112.0835 H21  178.6526
            X21 F21 1.0000 C21  90.0000 C22  180.0000
            0 1
        '''

        atom_name = [
          'F21', 'C1', 'C2', 'H11', 'H12', 'H13', 'H22', 'H23'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_b_carbon_zmatrix(self):

        zmatrix = '''\
                  H25  :1 DISTANCE :2  ANGLE  :3  DIHEDRAL
                  X21 H25 1.0000   :1  90.0000  :2   0.0000
                  C21 H25 1.1006 X21  90.0000   :1  180.0000
                  F21 C21 1.4174 H25 112.0835   :1  180.0000
                  C22 C21 1.5152 F21 109.5365 H25   100.0000
                  H21 C22 1.0996 C21 110.5586 F21   60.3866
                  H22 C22 1.0996 C21 110.5586 F21  -60.3866
                  H23 C22 1.1020 C21 109.4385 F21 -180.0000
                  H24 C21 1.1006 C22 112.0835 H21  -57.8794
                  X21 H25 1.0000 C21  90.0000 F21    0.0000
                  0 1
        '''

        atom_name = [
          'H22', 'C2', 'F21', 'C1' 'H11', 'H12', 'H13',  'H23'
        ]

        return textwrap.dedent(zmatrix), atom_name

