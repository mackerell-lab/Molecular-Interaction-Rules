#!/usr/bin/env python3
#
# Lennard-Jones-Drill-2: Dimethylamine
# ------------------------------------

# Imports
# -------

import textwrap

class Dimethylamine(object):

    def __init__(self):

        self.resi_name = 'DMAM'

    def get_monomer_a_species(self):

        '''

        Get the Monomer A Species

        '''

        monomer_a_species = {
            'N1': self.get_monomer_a_nitrogen_zmatrix(),
            'H1': self.get_monomer_a_nitrogen_hydrogen_zmatrix(),
            'H2': self.get_monomer_a_carbon_hydrogen_zmatrix()
        }

        return monomer_a_species

    def get_monomer_a_nitrogen_zmatrix(self):

        zmatrix = '''\
          N11
          C11  N11  1.4468
          H11  C11  1.0932  N11  113.8619
          H12  C11  1.0839  N11  109.6507  H11  121.4971
          H13  C11  1.0850  N11  109.4190  H11 -120.5584
          C12  N11  1.4468  C11  113.4553  H11   56.3756
          H14  C12  1.0850  N11  109.4190  C11   64.1828
          H15  C12  1.0839  N11  109.6507  C11 -177.8728
          H16  C12  1.0932  N11  113.8620  C11  -56.3756
          H17  N11  1.0005  C11  109.7680  H11  -66.8335
          X11  N11  1.0000  C11   90.0000  H11    0.0000
          0 1
        '''

        atom_name = [
            'N1', 'C1', 'H11', 'H12', 'H13', 'C2', 'H21', 'H22', 'H23', 'HN1'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_a_nitrogen_hydrogen_zmatrix(self):

        zmatrix = '''\
            H11
            X11 H11 1.0000
            N11 H11 1.0005 X11 90.0000
            C11 N11 1.4468 H11 113.0000 X11 30.0000
            H12 C11 1.0932 N11 113.8619 H11 0.8335
            H13 C11 1.0839 N11 109.6507 H11 121.4971
            H14 C11 1.0850 N11 109.4190 H11 -120.5584
            C12 N11 1.4468 H11 113.4553 X11 240.0
            H15 C12 1.0850 N11 109.4190 C11 64.1828
            H16 C12 1.0839 N11 109.6507 C11 -177.8728
            H17 C12 1.0932 N11 113.8620 C11 -56.3756
            0 1
        '''

        atom_name = [
          'HN1', 'N1', 'C1', 'H11', 'H12', 'H13', 'C2', 'H21', 'H22', 'H23',
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_a_carbon_hydrogen_zmatrix(self):

        zmatrix = '''\
          H11
          C11  H11  1.0932
          N11  C11   1.4468 H11  113.8619
          H12  C11  1.0839  N11  109.6507  H11  121.4971
          H13  C11  1.0850  N11  109.4190  H11 -120.5584
          C12  N11  1.4468  C11  113.4553  H11   56.3756
          H14  C12  1.0850  N11  109.4190  C11   64.1828
          H15  C12  1.0839  N11  109.6507  C11 -177.8728
          H16  C12  1.0932  N11  113.8620  C11  -56.3756
          H17  N11  1.0005  C11  109.7680  H11  -66.8335
          X11  H11  1.0000  C11   90.0000  N11    0.0000
          0 1
        '''

        atom_name = [
          'H11', 'C1', 'N1', 'H12', 'H13', 'C2', 'H21', 'H22', 'H23', 'HN1'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_b_species(self):

        monomer_b_species = {
          # 'H1': self.get_monomer_b_nitrogen_hydrogen_zmatrix(),
          'H2': self.get_monomer_b_carbon_hydrogen_zmatrix()
        }

        return monomer_b_species

    def get_monomer_b_nitrogen_hydrogen_zmatrix(self):

        zmatrix = '''\
           H21  :1   DISTANCE  :2 ANGLE  :3  DIHEDRAL
           X21 H21 1.0000  :1  90.0000   :2    0.0000
           N21 H21 1.0005 X21 90.0000    :1  180.0000
           C21 N21 1.4468 H21 113.0000 X21 30.0000
           H22 C21 1.0932 N21 113.8619 H21 0.8335
           H23 C21 1.0839 N21 109.6507 H21 121.4971
           H24 C21 1.0850 N21 109.4190 H21 -120.5584
           C22 N21 1.4468 H21 113.4553 X21 240.0
           H25 C22 1.0850 N21 109.4190 C21 64.1828
           H26 C22 1.0839 N21 109.6507 C21 -177.8728
           H27 C22 1.0932 N21 113.8620 C21 -56.3756
           0 1
        '''

        atom_names = [
          'HN1', 'N1', 'C1', 'H11', 'H12', 'H13', 'C2', 'H21', 'H22', 'H23',
        ]

        return textwrap.dedent(zmatrix), atom_names

    def get_monomer_b_carbon_hydrogen_zmatrix(self):

        zmatrix = '''\
          H21   :1   DISTANCE  :2  ANGLE     :3   DIHEDRAL
          X21  H21    1.0000   :1  90.0000     :2    0.0000
          C21  H21  1.0932     X21   90.0000   :2   180.0000
          N21  C21  1.4468 H21  113.8619       :1  180.0000
          H22  C21  1.0839  N21  109.6507  H21  121.4971
          H23  C21  1.0850  N21  109.4190  H21 -120.5584
          C22  N21  1.4468  C21  113.4553  H21   56.3756
          H24  C22  1.0850  N21  109.4190  C21   64.1828
          H25  C22  1.0839  N21  109.6507  C21 -177.8728
          H26  C22  1.0932  N21  113.8620  C21  -56.3756
          H27  N21  1.0005  C21  109.7680  H21  -66.8335
          0 1
        '''

        atom_names = [
          'H11', 'C1', 'N1', 'H12', 'H13', 'C2', 'H21', 'H22', 'H23', 'HN1'

        ]

        return textwrap.dedent(zmatrix), atom_names

