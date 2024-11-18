#!/usr/bin/env python3
#
# Lennard-Jones-Drill-2: EthylAmmonium
# ------------------------------------

# Imports
# -------

import textwrap

class EthylAmmonium(object):

    def __init__(self):

        self.resi_name = 'EAMM'

    def get_monomer_a_species(self):

        '''

        Get the Monomer A Species

        '''

        monomer_a_species = {
            'H1': self.get_monomer_a_nitrogen_hydrogen_zmatrix()
        }

        return monomer_a_species

    def get_monomer_a_nitrogen_hydrogen_zmatrix(self):

        zmatrix = '''\
            H11
            N11  H11   1.0870
            C11  N11   1.4558 H11  110.9666
            C12  C11   1.5226 N11  110.4162  H11  176.8461
            H12  C12  1.0870  C11  110.9666  N11   62.0187
            H13  C12  1.0870  C11  110.9666  N11  -57.6114
            H14  C12  1.0870  C11  110.9666  N11 -178.2677
            H15  C11  1.0851  C12  109.3856  H11  -63.1534
            H16  C11  1.0851  C12  109.3856  H11 -179.6728
            H17  N11  1.0870  C11  110.9666  C12  -64.9708
            H18  N11  1.0870  C11  -110.9666  C12  -120.0000
            X11  N11  1.0000  C11   90.0000  C12    0.0000
            1 1
        '''

        atom_name = [
            'HZ1', 'NZ', 'CE', 'C1', 'H11', 'H12', 'H13', 'HE1', 'HE2', 'HZ2', 'HZ3'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_b_species(self):

        monomer_b_species = {
        }

        return monomer_b_species

    def get_monomer_b_nitrogen_hydrogen_zmatrix(self):

      zmatrix = '''\
            H21   :1   DISTANCE     :2  ANGLE   :3   DIHEDRAL
            X21  H21    1.0000      :1  90.0000   :2    0.0000
            N21  H21   1.0870      X21   90.0000  :2  180.0000
            C21  N21   1.4558 H21  110.9666       :1  180.0000
            C22  C21   1.5226 N21  110.4162  H21  176.8461
            H22  C22  1.0870  C21  110.9666  N21   62.0187
            H23  C22  1.0870  C21  110.9666  N21  -57.6114
            H24  C22  1.0870  C21  110.9666  N21 -178.2677
            H25  C21  1.0851  C22  109.3856  H21  -63.1534
            H26  C21  1.0851  C22  109.3856  H21 -179.6728
            H27  N21  1.0870  C21  110.9666  C22  -64.9708
            H28  N21  1.0870  C21  -110.9666  C22  -120.0000
            X11  N11  1.0000  C11   90.0000  C12    0.0000
            1 1
        '''

      atom_name = [
        'HZ1', 'NZ', 'CE', 'C1', 'H11', 'H12', 'H13', 'HE1', 'HE2', 'HZ2', 'HZ3'
      ]

      return textwrap.dedent(zmatrix), atom_name



