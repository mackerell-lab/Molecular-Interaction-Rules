#!/usr/bin/env python3
#
# Lennard-Jones-Drill-2: EthylSulfanyl Phosphonic Acid
# -----------------------------------------------------

# Imports
# -------

import textwrap

class EthylSulfanylPhosphonicAcid(object):

  def __init__(self):

    self.resi_name = 'DMDS'

  def get_monomer_a_species(self):

    '''

    Get the Monomer A Species

    '''

    monomer_a_species = {
      'S1': self.get_monomer_a_sulphur_lone_pair_donor_zmatrix()
    }

    return monomer_a_species

  def get_monomer_b_species(self):

    '''

    Get the Monomer B Species

    '''

    monomer_b_species = {
      'S1': self.get_monomer_b_sulphur_lone_pair()
    }

    return monomer_b_species

  def get_monomer_a_sulphur_lone_pair_donor_zmatrix(self):

    zmatrix = '''\
            S11
            S12 S11  2.0847
            C11 S12  1.8292  S11  100.6848
            C12 S11  1.8292  S12  100.8205  C11   83.9593
            H11 C11  1.1017  S12  106.5343  S11  178.1127
            H12 C11  1.1017  S12  110.6848  S11   59.1606
            H13 C11  1.1017  S12  110.6848  S11  -63.3317
            H14 C12  1.1017  S11  110.6848  S12  -63.3317
            H15 C12  1.1017  S11  110.6848  S12   59.1606
            H16 C12  1.1017  S11  106.5343  S12  178.1127
            X11 S11  1.0000  C11   90.0000  C12  180.0000
            0 1
        '''

    atom_name = [
      'S2', 'S3', 'CM1', 'CM4', 'H1', 'H2', 'H3', 'H4', 'H5', 'H6'
    ]

    return textwrap.dedent(zmatrix), atom_name

  def get_monomer_b_sulphur_lone_pair(self):

    zmatrix = '''\
            S21  :1  DISTANCE  :2   90.0000  :3  180.0000
            S22 S21  2.0847     :1  126.2746  :2  180.0000
            C21 S22  1.8292  S21  100.6848   :1    0.0000
            C22 S21  1.8292  S22  100.8205  C21   83.9593
            H21 C21  1.1017  S22  106.5343  S21  178.1127
            H22 C21  1.1017  S22  110.6848  S21   59.1606
            H23 C21  1.1017  S22  110.6848  S21  -63.3317
            H24 C22  1.1017  S21  110.6848  S22  -63.3317
            H25 C22  1.1017  S21  110.6848  S22   59.1606
            H26 C22  1.1017  S21  106.5343  S22  178.1127
            X21 S21  1.0000  C21   90.0000  C22  180.0000
            0 1
        '''

    atom_name = [
      'S2', 'S3', 'CM1', 'CM4', 'H1', 'H2', 'H3', 'H4', 'H5', 'H6'
    ]

    return textwrap.dedent(zmatrix), atom_name

