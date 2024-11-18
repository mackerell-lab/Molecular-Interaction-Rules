#!/usr/bin/env python3
#
# Lennard-Jones-Drill-2: MethylPhosphate
# ---------------------------------------

# Imports
# -------

import textwrap

class MethylPhosphate(object):

    def __init__(self):

        self.resi_name = 'MP_0'

    def get_monomer_a_species(self):

        monomer_a_species = {
            'P1': self.get_phosphorous_zmatrix()
        }

        return monomer_a_species

    def get_phosphorous_zmatrix(self):

        zmatrix = '''\
          P11
          O11  P11  1.6437
          O12  P11  1.6437  O11   97.8486
          H11  O12  0.9705  P11  108.6226  O11  146.0940
          O13  P11  1.6286  O12  104.4218  H11 -106.7457
          C11  O13  1.4523  P11  117.8616  O11   51.1114
          H12  C11  1.0992  O13  110.2301  P11  -60.9713
          H13  C11  1.0992  O13  110.2301  P11   60.9714
          H14  C11  1.0957  O13  105.5916  P11 -179.9999
          O14  P11  1.5016  O12  117.4240  H11   19.5590
          H15  O11  0.9705  P11  108.6223  O12 -146.0925
          X11  P11  1.0000  O11   90.0000  O12  180.0000
          0 1
        '''

        atom_name = [
          'P1', 'O2', 'O3', 'H3', 'O1', 'C1', 'H11', 'H12', 'H13', 'O4', 'H2'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_b_species(self):

        monomer_b_species = {
            'P1': self.get_monomer_b_phosphorous_zmatrix()
        }

        return monomer_b_species

    def get_monomer_b_phosphorous_zmatrix(self):

        zmatrix = '''\
          P21   :1  DISTANCE   :2  ANGLE   :3    DIHEDRAL
          O21  P21  1.6437      :1  180.0000  :2   180.0000
          O22  P21  1.6437  O21   97.8486     :1     0.0000
          H21  O22  0.9705  P21  108.6226  O21  146.0940
          O23  P21  1.6286  O22  104.4218  H21 -106.7457
          C21  O23  1.4523  P21  117.8616  O21   51.1114
          H22  C21  1.0992  O23  110.2301  P21  -60.9713
          H23  C21  1.0992  O23  110.2301  P21   60.9714
          H24  C21  1.0957  O23  105.5916  P21 -179.9999
          O24  P21  1.5016  O22  117.4240  H21   19.5590
          H25  O21  0.9705  P21  108.6223  O22 -146.0925
          X21  P21  1.0000  O21   90.0000  O22  180.0000
          0 1
        '''


        atom_name = [
          'P1', 'O2', 'O3', 'H3', 'O1', 'C1', 'H11', 'H12', 'H13', 'O4', 'H2'
        ]

        return textwrap.dedent(zmatrix), atom_name

