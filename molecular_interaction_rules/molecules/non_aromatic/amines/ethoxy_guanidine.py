#!/usr/bin/env python3
#
# Lennard-Jones-Drill-2: Dimethylamine
# ------------------------------------

# Imports
# -------

import textwrap

class EthoxyGuanidine(object):

  def __init__(self):

    self.resi_name = 'SM086'

  def get_monomer_a_species(self):

      '''

      Get the Monomer A Species

      '''

      monomer_a_species = {
          'C1': self.get_monomer_a_carbon_zmatrix(),
          # 'H1': self.get_monomer_a_nitrogen__hydrogen_zmatrix()
      }

      return monomer_a_species

  def get_monomer_a_carbon_zmatrix(self):

      zmatrix = '''\
          C11
          N11   C11   1.3022
          O11   N11   1.4376    C11  108.4198
          N12   C11   1.3720    N11  126.2972  O11    6.4406
          H11   N12   0.9970    C11  114.2111  N11 -143.9845
          H12   N12   0.9964    C11  113.4884  N11  -11.5035
          N13   C11   1.3868    N11  119.2445  O11 -173.8444
          H13   N13   0.9998    C11  112.7731  N11 -131.4028
          H14   N13   0.9988    C11  111.5425  N11   -5.3372
          C12   O11   1.4008    N11  109.5895  C11  179.0476
          H15   C12   1.0852    O11  109.5408  N11   59.0744
          H16   C12   1.0852    O11  109.5624  N11  -59.0583
          C13   C12   1.5163    O11  107.7264  N11  179.9869
          H17   C13   1.0852    C12  110.6876  O11  -60.0019
          H18   C13   1.0852    C12  110.2510  O11 -179.9243
          H19   C13   1.0852    C12  110.7131  O11   60.1648
          X11   C11   1.0000    N11   90.0000  O11    0.0000
          0 1
      '''

      atom_name = [
        'CZ', 'NE', 'OD', 'NH1', 'HH11', 'HH12', 'NH2', 'HH21', 'HH22',
        'CG', 'HG1', 'HG2', 'CB', 'HB1', 'HB2', 'HB3'
      ]

      return textwrap.dedent(zmatrix), atom_name

  def get_monomer_a_nitrogen__hydrogen_zmatrix(self):


      zmatrix = '''\
          H11
          N12   H11   0.9964
          C11   N12   1.3720    H11  114.2111
          N11   C11   1.2644    N12  126.2972  H11 -143.9845
          O11   N11   1.3918    C11  110.5211  N12    6.4406
          H12   N12   0.9964    C11  113.4884  N11  -11.5035
          N13   C11   1.3868    N11  119.2445  O11 -173.8444
          H13   N13   0.9998    C11  112.7731  N11 -131.4028
          H14   N13   0.9988    C11  111.5425  N11   -5.3372
          C12   O11   1.4008    N11  109.5895  C11  179.0476
          H15   C12   1.0852    O11  109.5408  N11   59.0744
          H16   C12   1.0854    O11  109.5624  N11  -59.0583
          C13   C12   1.5163    O11  107.7264  N11  179.9869
          H17   C13   1.0845    C12  110.6876  O11  -60.0019
          H18   C13   1.0855    C12  110.2510  O11 -179.9243
          H19   C13   1.0846    C12  110.7131  O11   60.1648
          X11   H11   1.0000    N12   90.0000  C11    0.0000
          0 1
      '''

      atom_name = [
          'HH11', 'NH1', 'CZ', 'NE', 'OD', 'HH12', 'NH2', 'HH21', 'HH22',
          'CG', 'HG1', 'HG2', 'CB', 'HB1', 'HB2', 'HB3'
      ]

      return textwrap.dedent(zmatrix), atom_name

  def get_monomer_b_species(self):

      monomer_b_species = {
        'C1': self.get_monomer_b_carbon_zmatrix(),
      }

      return monomer_b_species

  def get_monomer_b_carbon_zmatrix(self):

    zmatrix = '''\
          C21    :1   DISTANCE   :2  ANGLE       :3   DIHEDRAL
          N21   C21   1.2644     :1  236.0000    :2     0.0000
          O21   N21   1.3918    C21  110.5211    :1    0.0000
          N22   C21   1.3720    N21  126.2972  O21    6.4406
          H21   N22   0.9970    C21  114.2111  N21 -143.9845
          H22   N22   0.9964    C21  113.4884  N21  -11.5035
          N23   C21   1.3868    N21  119.2445  O21 -173.8444
          H23   N23   0.9998    C21  112.7731  N21 -131.4028
          H24   N23   0.9988    C21  111.5425  N21   -5.3372
          C22   O21   1.4008    N21  109.5895  C21  179.0476
          H25   C22   1.0852    O21  109.5408  N21   59.0744
          H26   C22   1.0852    O21  109.5624  N21  -59.0583
          C23   C22   1.5163    O21  107.7264  N21  179.9869
          H27   C23   1.0852    C22  110.6876  O21  -60.0019
          H28   C23   1.0852    C22  110.2510  O21 -179.9243
          H29   C23   1.0852    C22  110.7131  O21   60.1648
          0 1
      '''

    atom_name = [
      'CZ', 'NE', 'OD', 'NH1', 'HH11', 'HH12', 'NH2', 'HH21', 'HH22',
      'CG', 'HG1', 'HG2', 'CB', 'HB1', 'HB2', 'HB3'
    ]

    return textwrap.dedent(zmatrix), atom_name
