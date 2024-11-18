#!/usr/bin/env python3
#
# Lennard-Jones-Drill-2: Phenol
# -----------------------------

# Imports
# -------
import textwrap

class Phenol(object):

  def __init__(self):

      self.resi_name = 'phen'

  def get_monomer_a_species(self):


    monomer_a_species = {
        'RC1': self.monomer_a_aromatic_zmatrix(),
        'O1': self.monomer_a_oxygen_zmatrix()
    }

    return monomer_a_species

  def get_monomer_b_species(self):

      monomer_b_species = {
        'RC1': self.monomer_b_aromatic_zmatrix(),
        'O1': self.monomer_b_oxygen_zmatrix()
      }

      return monomer_b_species

  def monomer_a_aromatic_zmatrix(self):

      zmatrix = '''\
          X11
          C11  X11  1.3940
          C12  C11  1.4080 X11   60.0000
          C13  C12  1.4080 C11  120.0000 X11    0.0000
          C14  C13  1.4080 C12  120.0000 C11    0.0000
          C15  C14  1.4080 C13  120.0000 C12    0.0000
          C16  C15  1.4080 C14  120.0000 C13    0.0000
          O11  C11  1.3820 C12  240.0000 C13    0.0000
          H11  O11  0.9560 C11  109.0000 C12  180.0000
          H12  C12  1.0939 C11  120.0000 C13  180.0000
          H13  C13  1.0939 C12  120.0000 C11  180.0000
          H14  C14  1.0939 C13  120.0000 C12  180.0000
          H15  C15  1.0939 C14  120.0000 C13  180.0000
          H16  C16  1.0939 C15  120.0000 C11  180.0000
          0 1
      '''

      atom_name = [
        'HH', 'CZ', 'CE1', 'CD1', 'CG', 'CD2', 'CE2', 'OH', 'HE1', 'HD1', 'HG', 'HD2', 'HE2'
      ]

      return textwrap.dedent(zmatrix), atom_name

  def monomer_a_oxygen_zmatrix(self):

    zmatrix = '''\
        O11
        H11  O11  0.9560
        X11  H11  1.0000 O11   90.000
        C11  O11  1.3820 H11  109.0000  X11   0.0000
        C12  C11  1.4080 O11   120.0000 H11   0.0000
        C13  C12  1.4080 C11   120.0000 O11  180.0000
        C14  C13  1.4080 C12  120.0000 C11    0.0000
        C15  C14  1.4080 C13  120.0000 C12    0.0000
        C16  C15  1.4080 C14  120.0000 C13    0.0000
        H12  C12  1.0939 C11  120.0000 C13  180.0000
        H13  C13  1.0939 C12  120.0000 C11  180.0000
        H14  C14  1.0939 C13  120.0000 C12  180.0000
        H15  C15  1.0939 C14  120.0000 C13  180.0000
        H16  C16  1.0939 C15  120.0000 C11  180.0000
        0 1
      '''

    atom_name = [
      'OH', 'CZ', 'CE1', 'CD1', 'CG', 'CD2', 'CE2', 'HH', 'HE1', 'HD1', 'HG', 'HD2', 'HE2'
    ]

    return textwrap.dedent(zmatrix), atom_name

  def monomer_b_aromatic_zmatrix(self):

      zmatrix = '''\
            X21   :1  DISTANCE  :2  ANGLE  :3   90.0000
            C21  X21  1.3940  :1  90.0000   :2  180.0000
            C22  C21  1.4080 X21   60.0000  :1   90.0000
            C23  C22  1.4080 C21  120.0000 X21    DIHEDRAL
            C24  C23  1.4080 C22  120.0000 C21    0.0000
            C25  C24  1.4080 C23  120.0000 C22    0.0000
            C26  C25  1.4080 C24  120.0000 C23    0.0000
            O21  C21  1.3820 C22  240.0000 C23    0.0000
            H21  O21  0.9560 C21  109.0000 C22  180.0000
            H22  C22  1.0939 C21  120.0000 C23  180.0000
            H23  C23  1.0939 C22  120.0000 C21  180.0000
            H24  C24  1.0939 C23  120.0000 C22  180.0000
            H25  C25  1.0939 C24  120.0000 C23  180.0000
            H26  C26  1.0939 C25  120.0000 C21  180.0000
            0 1
        '''

      atom_name = [
        'CZ', 'CE1', 'CD1', 'CG', 'CD2', 'CE2', 'OH', 'HH', 'HE1', 'HD1', 'HG', 'HD2', 'HE2'
      ]

      return textwrap.dedent(zmatrix), atom_name

  def monomer_b_oxygen_zmatrix(self):

    zmatrix = '''\
        O21   :1  DISTANCE  :2  ANGLE  :3  180.0000
        H21  O21  0.9560   :1  115.2746    :2  DIHEDRAL
        X21  H21  1.0000  O21   90.0000    :1   180.0000
        C21  O21  1.3820  H21  115.2746  X21    0.0000
        C22  C21  1.4080 O21  -240.0000  H21    0.0000
        C23  C22  1.4080 C21  -120.0000 O21    0.0000
        C24  C23  1.4080 C22  120.0000 C21    0.0000
        C25  C24  1.4080 C23  120.0000 C22    0.0000
        C26  C25  1.4080 C24  120.0000 C23    0.0000
        H22  C22  1.0939 C21  120.0000 C23  180.0000
        H23  C23  1.0939 C22  120.0000 C21  180.0000
        H24  C24  1.0939 C23  120.0000 C22  180.0000
        H25  C25  1.0939 C24  120.0000 C23  180.0000
        H26  C26  1.0939 C25  120.0000 C21  180.0000
        X21  O21  1.0000 C21   90.0000 C22    0.0000
        0 1
      '''

    atom_name = [
      'OH', 'CZ', 'CE1', 'CD1', 'CG', 'CD2', 'CE2', 'HH', 'HE1', 'HD1', 'HG', 'HD2', 'HE2'
    ]

    return textwrap.dedent(zmatrix), atom_name

