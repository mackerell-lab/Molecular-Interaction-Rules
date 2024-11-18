#!/usr/bin/env python3
#
# Lennard-Jones-Drill-2: Four Pyridinone
# --------------------------------------

# Imports
# -------
import textwrap

class FourPyridinone(object):

  def __init__(self):

    self.resi_name = '4PYO'

  def get_monomer_a_species(self):

    monomer_a_species = {
      'O1': self.monomer_a_oxygen_zmatrix(),
      'RC1': self.monomer_a_aromatic_zmatrix()
    }

    return monomer_a_species

  def get_monomer_b_species(self):

    monomer_b_species = {
      'O1': self.monomer_b_oxygen_zmatrix(),
      'RC1': self.monomer_b_aromatic_zmatrix()
    }

    return monomer_b_species

  def monomer_a_oxygen_zmatrix(self):

      zmatrix = '''\
          O11
          C11  O11  1.2500
          C12  C11  1.4671   O11  123.1498
          C13  C12  1.3741   C11  121.7979 O11 -180.0000
          N11  C13  1.3754   C12  120.9490 C11   -0.0000
          C14  N11  1.3754   C13  120.8056 C12   -0.0000
          C15  C14  1.3741   N11  120.9490 C13   -0.0000
          H11  C13  1.0920   C12  123.1203 C11 -180.0000
          H12  N11  1.0123   C13  119.5971 C12 -180.0000
          H13  C14  1.0920   N11  115.9305 C13 -180.0000
          H14  C15  1.0924   C11  118.6385 O11   -0.0000
          H15  C12  1.0924   C11  118.6385 O11   -0.0000
          X11  O11  1.0000   C11   90.0000 C12    0.0000
          0 1
      '''

      atom_name = [
        'O1', 'C1', 'C22', 'C32', 'N1', 'C31', 'C21', 'H32', 'H1', 'H31' ,'H21', 'H22',
      ]

      return textwrap.dedent(zmatrix), atom_name

  def monomer_a_aromatic_zmatrix(self):

      zmatrix = '''\
          X11
          C12  X11 1.2940
          C11  C12  1.4671   X11  60.0000
          O11  C11  1.2500   C12  123.1498 X11  180.0000
          C13  C12  1.3741   C11  121.7979 O11 -180.0000
          N11  C13  1.3754   C12  120.9490 C11   -0.0000
          C14  N11  1.3754   C13  120.8056 C12   -0.0000
          C15  C14  1.3741   N11  120.9490 C13   -0.0000
          H11  C13  1.0920   C12  123.1203 C11 -180.0000
          H12  N11  1.0123   C13  119.5971 C12 -180.0000
          H13  C14  1.0920   N11  115.9305 C13 -180.0000
          H14  C15  1.0924   C11  118.6385 O11   -0.0000
          H15  C12  1.0924   C11  118.6385 O11   -0.0000
          0 1
      '''

      atom_name = [
        'C22', 'C1', 'O1', 'C32', 'N1', 'C31', 'C21', 'H32', 'H1', 'H31' ,'H21', 'H22',
      ]

      return textwrap.dedent(zmatrix), atom_name

  def monomer_b_oxygen_zmatrix(self):

      zmatrix = '''\
          O21   :1  DISTANCE   :2   ANGLE      :3   90.0000
          X21  O21  1.0000     :1   90.0000    :2    0.0000
          C21  O21  1.2500     X21  90.0000    :1  DIHEDRAL
          C22  C21  1.3986     O21  120.0000   :1    0.0000
          C23  C22  1.3741   C21  121.7979 O21 -180.0000
          N21  C23  1.3754   C22  120.9490 C21   -0.0000
          C24  N21  1.3754   C23  120.8056 C22   -0.0000
          C25  C24  1.3741   N21  120.9490 C23   -0.0000
          H21  C23  1.0920   C22  123.1203 C21 -180.0000
          H22  N21  1.0123   C23  119.5971 C22 -180.0000
          H23  C24  1.0920   N21  115.9305 C23 -180.0000
          H24  C25  1.0924   C21  118.6385 O21   -0.0000
          H25  C22  1.0924   C21  118.6385 O21   -0.0000
          X21  O21  1.0000   C21   90.0000 C22    0.0000
          0 1
      '''

      atom_name = [
        'O1', 'C1', 'C22', 'C32', 'N1', 'C31', 'C21', 'H32', 'H1', 'H31' ,'H21', 'H22',
      ]

      return textwrap.dedent(zmatrix), atom_name

  def monomer_b_aromatic_zmatrix(self):

      zmatrix = '''\
          X21   :1  DISTANCE  :2   ANGLE   :3     90.0000
          C22  X21 1.2940     :1   90.0000    :2   0.0000
          C21  C22  1.4671   X21  60.0000     :1  DIHEDRAL
          O21  C21  1.2500   C22  123.1498 X21    0.0000
          C23  C22  1.3741   C21  121.7979 O21 -180.0000
          N21  C23  1.3754   C22  120.9490 C21   -0.0000
          C24  N21  1.3754   C23  120.8056 C22   -0.0000
          C25  C24  1.3741   N21  120.9490 C23   -0.0000
          H21  C23  1.0920   C22  123.1203 C21 -180.0000
          H22  N21  1.0123   C23  119.5971 C22 -180.0000
          H23  C24  1.0920   N21  115.9305 C23 -180.0000
          H24  C25  1.0924   C21  118.6385 O21   -0.0000
          H25  C22  1.0924   C21  118.6385 O21   -0.0000
          0 1
      '''

      atom_name = [
        'C22', 'C1', 'O1', 'C32', 'N1', 'C31', 'C21', 'H32', 'H1', 'H31' ,'H21', 'H22',
      ]

      return textwrap.dedent(zmatrix), atom_name

