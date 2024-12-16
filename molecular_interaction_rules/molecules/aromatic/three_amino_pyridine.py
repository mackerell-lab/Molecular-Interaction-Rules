#!/usr/bin/env python3
#
# Lennard-Jones-Drill-2: Pyridine
# -------------------------------

# Imports
# -------
import textwrap

class ThreeAminoPyridine(object):

  def __init__(self):

    self.resi_name = '3APY'

  def get_monomer_a_species(self):

    '''

    Get the Monomer A Species

    '''

    monomer_a_species = {
        'H1': self.monomer_a_nitrogen_zmatrix()
    }

    return monomer_a_species

  def get_monomer_b_species(self):

    monomer_b_species = {
      'H1': self.monomer_b_nitrogen_zmatrix()
    }

    return monomer_b_species

  def monomer_a_nitrogen_zmatrix(self):

      zmatrix = '''\
          H11
          N11  H11  1.0162
          C11  N11  1.4056   H11  114.1487
          C12  C11  1.4145   N11  120.8173   H11 -153.4115
          N12  C12  1.3483   C11  124.2902   N11 -176.4540
          C13  N12  1.3530   C12  117.3356   C11   -0.0680
          C14  C13  1.4059   N12  122.9426   C12   -0.0449
          C15  C14  1.4029   C13  119.2363   N12    0.1804
          H11  C13  1.0940   N12  116.2104   C12 -179.9682
          H12  C14  1.0936   C13  120.2082   N12 -179.8993
          H13  C15  1.0952   C11  120.3578   N11   -3.1214
          H14  C12  1.0972   C11  119.5461   N11    3.3597
          H15  N11  1.0164   C11  114.0036   C12  -24.8625
          X11  H11  1.0000   N11   90.0000   C11  180.0000
          0 1
        '''

      atom_name = [
          'HD12', 'ND1', 'CD1', 'CE1', 'NZ', 'CE2', 'CD2', 'CG', 'HE2', 'HD2', 'HG', 'HE1', 'HD11'
      ]

      return textwrap.dedent(zmatrix), atom_name

  def monomer_b_nitrogen_zmatrix(self):

      zmatrix = '''\
            H21  X11  DISTANCE   :2   ANGLE    :3   0.0000
            X21   H21  1.0000    :1   90.0000   :2   DIHEDRAL
            N21  H21  1.0162    X21   90.0000  :1  180.0000
            C21  N21  1.4056   H21  114.1487   X21   0.0000
            C22  C21  1.4145   N21  120.8173   H21 -153.4115
            N22  C22  1.3483   C21  124.2902   N21 -176.4540
            C23  N22  1.3530   C22  117.3356   C21   -0.0680
            C24  C23  1.4059   N22  122.9426   C22   -0.0449
            C25  C24  1.4029   C23  119.2363   N22    0.1804
            H21  C23  1.0940   N22  116.2104   C22 -179.9682
            H22  C24  1.0936   C23  120.2082   N22 -179.8993
            H23  C25  1.0952   C21  120.3578   N21   -3.1214
            H24  C22  1.0972   C21  119.5461   N21    3.3597
            H25  N21  1.0164   C21  114.0036   C22  -24.8625
            X21  H21  1.0000   N21   90.0000   C21  180.0000
            0 1
      '''

      atom_name = [
        'HD12', 'ND1', 'CD1', 'CE1', 'NZ', 'CE2', 'CD2', 'CG', 'HE2', 'HD2', 'HG', 'HE1', 'HD11'
      ]

      return textwrap.dedent(zmatrix), atom_name

