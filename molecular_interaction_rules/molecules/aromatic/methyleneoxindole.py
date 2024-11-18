#!/usr/bin/env python3
#
# Lennard-Jones-Drill-2: Methyleneoxindole
# ----------------------------------------

# Imports
# -------
import textwrap

class Methyleneoxindole(object):

  def __init__(self):

      self.resi_name = 'MEOI'

  def get_monomer_a_species(self):

    monomer_a_species = {
      'RC1': self.monomer_a_aromatic_zmatrix()
    }

    return monomer_a_species

  def get_monomer_b_species(self):

    monomer_b_species = {
        'RC1': self.monomer_b_aromatic_zmatrix()
    }

    return monomer_b_species

  def monomer_a_aromatic_zmatrix(self):

      zmatrix = '''\
          X11
          C11        X11    1.0122
          C12        C11    1.5142      X11   59.0000
          O11        C11    1.2283      C12  128.7573      X11    180.0000
          C13        C12    1.3546      C11  122.3840      O11   -0.0000
          H11        C13    1.0929      C12  121.6798      C11 -180.0000
          H12        C13    1.0933      C12  119.3351      C11   -0.0000
          C14        C12    1.4699      C11  106.6435      O11 -180.0000
          C15        C14    1.4034      C12  132.8041      C11 -180.0000
          H13        C15    1.0942      C14  120.6716      C12   -0.0000
          C16        C15    1.4088      C14  118.6956      C12 -180.0000
          H14        C16    1.0933      C15  119.8575      C14 -180.0000
          C17        C16    1.4111      C15  120.5682      C14   -0.0000
          H15        C17    1.0938     C16  119.5417      C15 -180.0000
          C18        C17    1.4107     C16  121.3126      C15   -0.0000
          H16        C18    1.0936     C17  120.9423     C16 -180.0000
          C19        C14    1.4189      C12  107.1110      C13 -180.0000
          N11        C11    1.3946      C12  105.1443    C13 -180.0000
          H17        N11    1.0139      C11  122.2020      O11   -0.0000
          0 1
      '''

      atom_name = [
        'NE1', 'CE2', 'CZ2', 'CH2', 'CZ3', 'CE3', 'CD2', 'CG', 'CD1',
        'HE1', 'HZ2', 'HH2', 'HZ3', 'HE3', 'HG', 'HD1'
      ]

      return textwrap.dedent(zmatrix), atom_name

  def monomer_b_aromatic_zmatrix(self):

    zmatrix = '''\
          X21        :1  DISTANCE       :2   ANGLE       :3   90.0000
          C21        X21    1.0122      :1   90.0000     :2  180.0000
          C22        C21    1.5142      X21   59.0000    :1   90.0000
          O21        C21    1.2283      C22  128.7573      X21    DIHEDRAL
          C23        C22    1.3546      C21  122.3840      O21   -0.0000
          H21        C23    1.0929      C22  121.6798      C21 -180.0000
          H22        C23    1.0933      C22  119.3351      C21   -0.0000
          C24        C22    1.4699      C21  106.6435      O21 -180.0000
          C25        C24    1.4034      C22  132.8041      C21 -180.0000
          H23        C25    1.0942      C24  120.6716      C22   -0.0000
          C26        C25    1.4088      C24  118.6956      C22 -180.0000
          H24        C26    1.0933      C25  119.8575      C24 -180.0000
          C27        C26    1.4111      C25  120.5682      C24   -0.0000
          H25        C27    1.0938     C26  119.5417      C25 -180.0000
          C28        C27    1.4107     C26  121.3126      C25   -0.0000
          H26        C28    1.0936     C27  120.9423     C26 -180.0000
          C29        C24    1.4189      C22  107.1110      C23 -180.0000
          N21        C21    1.3946      C22  105.1443    C23 -180.0000
          H27        N21    1.0139      C21  122.2020      O21   -0.0000
          0 1
    '''

    atom_name = [
        'NE1', 'CE2', 'CZ2', 'CH2', 'CZ3', 'CE3', 'CD2', 'CG', 'CD1',
        'HE1', 'HZ2', 'HH2', 'HZ3', 'HE3', 'HG', 'HD1'
    ]

    return textwrap.dedent(zmatrix), atom_name



