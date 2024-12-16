#!/usr/bin/env python3
#
# Lennard-Jones-Drill-2: Uracil
# -----------------------------

# Imports
# -------
import textwrap

class Uracil(object):

  def __init__(self):

    self.resi_name = 'URAC'

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
          C11 O11  1.2314
          C12 C11  1.4628  O11  126.1129
          C13 C12  1.3634  C11  119.6473  O11 -180.0000
          N11 C13  1.3808  C12  121.7593  C11   -0.0000
          C14 N11  1.3939  C13  123.6682  C12   -0.0000
          O12 C14  1.2279  N11  122.9749  C13 -180.0000
          N12 C11  1.4125  C12  113.5676  C13   -0.0000
          H11 C13  1.0916  C12  122.6427  C11 -180.0000
          H12 N11  1.0139  C13  121.2538  C12 -180.0000
          H13 N12  1.0183  C11  116.2101  O11   -0.0000
          H14 C12  1.0888  C11  118.8137  O11   -0.0000
          X11 O11  1.0000  C11   90.0000  C12  180.0000
          0 1
      '''

      atom_name = [
         'O4', 'C4' , 'C5', 'C6', 'N1', 'C2', 'O2', 'N3', 'H6', 'H1', 'H3', 'H5'
      ]

      return textwrap.dedent(zmatrix), atom_name

  def monomer_a_aromatic_zmatrix(self):

      zmatrix = '''\
          X11
          C12 X11  1.2940
          C11 C12  1.4628  X11   59.0000
          O11 C11  1.2940  C12  126.1129  X11  180.0000
          C13 C12  1.3634  C11  119.6473  O11 -180.0000
          N11 C13  1.3808  C12  121.7593  C11   -0.0000
          C14 N11  1.3939  C13  123.6682  C12   -0.0000
          O12 C14  1.2279  N11  122.9749  C13 -180.0000
          N12 C11  1.4125  C12  113.5676  C13   -0.0000
          H11 C13  1.0916  C12  122.6427  C11 -180.0000
          H12 N11  1.0139  C13  121.2538  C12 -180.0000
          H13 N12  1.0183  C11  116.2101  O11   -0.0000
          H14 C12  1.0888  C11  118.8137  O11   -0.0000
          X11 O11  1.0000  C11   90.0000  C12  180.0000
          0 1
      '''

      atom_name = [
        'C5', 'C4', 'O4','C6', 'N1', 'C2', 'O2', 'N3', 'H6', 'H1', 'H3', 'H5'
      ]

      return textwrap.dedent(zmatrix), atom_name

  def monomer_b_oxygen_zmatrix(self):

      zmatrix = '''\
            O21  :1  DISTANCE  :2    ANGLE   :3   90.0000
            C21 O21  1.2314    :1  180.0000   :2  DIHEDRAL
            C22 C21  1.4628  O21  -126.1129   :1    0.0000
            C23 C22  1.3634  C21  119.6473  O21 -180.0000
            N21 C23  1.3808  C22  121.7593  C21   -0.0000
            C24 N21  1.3939  C23  123.6682  C22   -0.0000
            O22 C24  1.2279  N21  122.9749  C23 -180.0000
            N22 C21  1.4125  C22  113.5676  C23   -0.0000
            H21 C23  1.0916  C22  122.6427  C21 -180.0000
            H22 N21  1.0139  C23  121.2538  C22 -180.0000
            H23 N22  1.0183  C21  116.2101  O21   -0.0000
            H24 C22  1.0888  C21  118.8137  O21   -0.0000
            X21 O21  1.0000  C21   90.0000  C22  180.0000
            0 1
        '''

      atom_name = [
        'O4', 'C4' , 'C5', 'C6', 'N1', 'C2', 'O2', 'N3', 'H6', 'H1', 'H3', 'H5'
      ]

      return textwrap.dedent(zmatrix), atom_name

  def monomer_b_aromatic_zmatrix(self):

      zmatrix = '''\
            X21  :1  DISTANCE  :2   ANGLE  :3   90.0000
            C22 X21  1.2940    :1   90.0000   :2  180.0000
            C21 C22  1.4628  X21   59.0000    :1   90.0000
            O21 C21  1.2940  C22  126.1129  X21  DIHEDRAL
            C23 C22  1.3634  C21  119.6473  O21 -180.0000
            N21 C23  1.3808  C22  121.7593  C21   -0.0000
            C24 N21  1.3939  C23  123.6682  C22   -0.0000
            O22 C24  1.2279  N21  122.9749  C23 -180.0000
            N22 C21  1.4125  C22  113.5676  C23   -0.0000
            H21 C23  1.0916  C22  122.6427  C21 -180.0000
            H22 N21  1.0139  C23  121.2538  C22 -180.0000
            H23 N22  1.0183  C21  116.2101  O21   -0.0000
            H24 C22  1.0888  C21  118.8137  O21   -0.0000
            0 1
      '''

      atom_name = [
        'C5', 'C4', 'O4','C6', 'N1', 'C2', 'O2', 'N3', 'H6', 'H1', 'H3', 'H5'
      ]

      return textwrap.dedent(zmatrix), atom_name

