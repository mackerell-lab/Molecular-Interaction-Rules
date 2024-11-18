#!/usr/bin/env python3
#
# Lennard-Jones-Drill-2: Tribromoethane
# -------------------------------------

# Imports
# -------

import textwrap

class Tribromoethane(object):

  def __init__(self):

    self.resi_name = 'TBRE'

  def get_monomer_a_species(self):

    '''

    Get the Monomer A Species

    '''

    monomer_a_species = {
        'BR1': self.get_monomer_a_bromo_zmatrix()
    }

    return monomer_a_species

  def get_monomer_a_bromo_zmatrix(self):

      zmatrix = '''\
          BR11
          C11   BR11   1.9400
          C12   C11    1.5054   BR11  110.1466
          H11   C12    1.0996   C11  109.4854    BR11  -58.6258
          H12   C12    1.0996   C11  109.4854    BR11 -178.6928
          H13   C12    1.0996   C11  109.4854    BR11   61.4412
          BR12   C11    1.9400   C12  115.1421    H11 -180.0000
          BR13   C11   1.9400   C12  115.1421    H11   58.6258
          X11   BR11   1.0000   C11   90.0000    C12  180.0000
          0 1
       '''

      atom_name = [
        'BR12', 'C2', 'C1', 'H11', 'H12', 'H13', 'BR11', 'BR13'
      ]

      return textwrap.dedent(zmatrix), atom_name

  def get_monomer_b_species(self):

    monomer_b_species = {
      'BR1': self.get_monomer_b_bromo_zmatrix(),
    }

    return monomer_b_species

  def get_monomer_b_bromo_zmatrix(self):

    zmatrix = '''\
      BR21   :1    DISTANCE    :2 ANGLE  :3  DIHEDRAL
      X21   BR21   1.0000      :1  90.0000     :2   0.0000
      C21   BR21   1.9400    X21  90.0000    :1  180.0000
      C22   C21    1.5054   BR21  110.1466   :1  180.0000
      H21   C22    1.0996   C21  109.4854    BR21  -58.6258
      H22   C22    1.0996   C21  109.4854    BR21 -178.6928
      H23   C22    1.0996   C21  109.4854    BR21   61.4412
      BR22   C21   1.9400   C22  115.1421    H21 -180.0000
      BR23   C21   1.9400   C22  115.1421    H21   58.6258
      0 1
    '''

    atom_name = [
      'BR12', 'C2', 'C1', 'H11', 'H12', 'H13', 'BR11', 'BR13'
    ]

    return zmatrix, atom_name



