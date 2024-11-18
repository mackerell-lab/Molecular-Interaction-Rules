#!/usr/bin/env python3
#
# Lennard-Jones-Drill-2: Amidinium
# --------------------------------

# Imports
# -------

import textwrap

class Amidinium(object):

  def __init__(self):

    self.resi_name = 'AMDN'

  def get_monomer_a_species(self):

    monomer_a_species = {
      'H1': self.get_monomer_a_nitrogen_zmatrix(),
    }

    return monomer_a_species

  def get_monomer_a_nitrogen_zmatrix(self):

        zmatrix = '''\
          H16
          N11  H16  1.0158
          C11  N11  1.3802   H16  113.5777
          N12  C11  1.2582   N11  119.8884  H16   13.7994
          H11  N12  1.0019   C11  111.4728  N11  175.6979
          H12  N12  1.0019   C11  111.4728  N11   60.0000
          C12  C11  1.5109   N12  126.4061  H11   -1.4460
          H12  C12  1.0980   C11  110.6754  N11  -57.1481
          H13  C12  1.0980   C11  110.6754  N11   61.7100
          H14  C12  1.0980   C11  110.6754  N11 -177.9849
          H15  N11  0.9961   C11  116.9658  N12  149.2893
          X11  H16  1.0000   N11   90.0000  C11    0.0000
          1 1
        '''

        atom_name = [
          'H11', 'N6', 'C1', 'N9', 'H7', 'H8', 'C2', 'H3', 'H4', 'H5', 'H10',
        ]

        return textwrap.dedent(zmatrix), atom_name

  def get_monomer_b_species(self):

        '''

        Get the Monomer B Species

        '''

        monomer_b_species = {
          'H1': self.get_monomer_b_nitrogen_zmatrix(),
        }

        return monomer_b_species

  def get_monomer_b_nitrogen_zmatrix(self):

        zmatrix = '''\
                H26  :1   DISTANCE  :2  ANGLE      :3  DIHEDRAL
                X21  H26  1.0000    :1  90.0000     :2    0.0000
                N21  H26  1.0158   X21  90.0000     :2  180.0000
                C21  N21  1.3802   H26  113.5777    :1  180.0000
                N22  C21  1.2582   N21  119.7147  H26   13.7994
                H21  N22  1.0019   C21  111.4728  N21  175.6979
                H22  N22  1.0019   C21  111.4728  N21   60.0000
                C22  C21  1.5109   N22  126.4061  H21   -1.4460
                H23  C22  1.0847   C21  110.6754  N21  -57.1481
                H24  C22  1.0863   C21  110.6754  N21   61.7100
                H25  C22  1.0818   C21  110.6754  N21 -177.9849
                H27  N21  0.9961   C21  116.9658  N22  149.2893
                1 1
            '''

        atom_name = [
            'H11', 'N6', 'C1', 'N9', 'H7', 'H8', 'C2', 'H3', 'H4', 'H5', 'H10',
        ]

        return textwrap.dedent(zmatrix), atom_name

