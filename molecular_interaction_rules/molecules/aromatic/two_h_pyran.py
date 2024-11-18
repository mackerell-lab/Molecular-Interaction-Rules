#!/usr/bin/env python3
#
# Lennard-Jones-Drill-2: Pyran
# ----------------------------

# Imports
# -------
import textwrap

class TwoHPyran(object):

  def __init__(self):

    self.resi_name = 'PY02'

  def get_monomer_a_species(self):

    monomer_a_species = {
        'O1': self.get_monomer_a_oxygen_zmatrix(),
        'RC1': self.monomer_a_aromatic_zmatrix()
    }

    return monomer_a_species

  def get_monomer_b_species(self):

    monomer_b_species = {
      'O1': self.get_monomer_b_oxygen_zmatrix(),
      'RC1': self.monomer_b_aromatic_zmatrix()
    }

    return monomer_b_species

  def get_monomer_a_oxygen_zmatrix(self):

    zmatrix = '''\
        O11
        C11   O11   1.3662
        C12   C11   1.3665   O11  122.6854
        C13   C12   1.4610   C11  117.7425   O11    8.3001
        C14   C13   1.3618   C12  118.3221   C11  -19.0238
        C15   C14   1.5059   C13  116.1256   C12   -7.3457
        H11   C13   1.0938   C12  119.8205   C11  155.4737
        H12   C14   1.0931   C13  122.9298   C12  173.7068
        H13   C15   1.0974   O11  105.0592   C11 -174.6970
        H14   C15   1.1101   O11  106.9039   C11   68.8999
        H15   C12   1.0910   C11  119.1790   O11  178.6725
        H16   C11   1.0927   C12  124.3425   C13 -165.4041
        X11   O11   1.0000   C11   90.0000   C12   180.0000
        0 1
    '''

    atom_name = [
       'O1', 'C6', 'C5', 'C4', 'C3', 'C2', 'H41', 'H31', 'H21', 'H22', 'H51', 'H61'
    ]

    return textwrap.dedent(zmatrix), atom_name

  def monomer_a_aromatic_zmatrix(self):

    zmatrix = '''\
        X11
        O11   X11   1.4950
        C11   O11   1.3662   X11   59.0000
        C12   C11   1.3665   O11  122.6854   X11    0.0000
        C13   C12   1.4610   C11  117.7425   O11    8.3001
        C14   C13   1.3618   C12  118.3221   C11  -19.0238
        C15   C14   1.5059   C13  116.1256   C12   -7.3457
        H11   C13   1.0938   C12  119.8205   C11  155.4737
        H12   C14   1.0931   C13  122.9298   C12  173.7068
        H13   C15   1.0974   O11  105.0592   C11 -174.6970
        H14   C15   1.1101   O11  106.9039   C11   68.8999
        H15   C12   1.0910   C11  119.1790   O11  178.6725
        H16   C11   1.0927   C12  124.3425   C13 -165.4041
        0 1
    '''

    atom_name = [
      'O1', 'C6', 'C5', 'C4', 'C3', 'C2', 'H41', 'H31', 'H21', 'H22', 'H51', 'H61'
    ]

    return textwrap.dedent(zmatrix), atom_name

  def get_monomer_b_oxygen_zmatrix(self):

    zmatrix = '''\
        O21    :1  DISTANCE    :2  ANGLE      :3   90.0000
        C21   O21   1.3662     :1  236.0000   :2   DIHEDRAL
        C22   C21   1.3665   O21  -122.6854   :1    0.0000
        C23   C22   1.4610   C21  117.7425   O21    8.3001
        C24   C23   1.3618   C22  118.3221   C21  -19.0238
        C25   C24   1.5059   C23  116.1256   C22   -7.3457
        H21   C23   1.0938   C22  119.8205   C21  155.4737
        H22   C24   1.0931   C23  122.9298   C22  173.7068
        H23   C25   1.0974   O21  105.0592   C21 -174.6970
        H24   C25   1.1101   O21  106.9039   C21   68.8999
        H25   C22   1.0910   C21  119.1790   O21  178.6725
        H26   C21   1.0927   C22  124.3425   C23 -165.4041
        X21   O21   1.0000   C21   90.0000   C22    0.0000
        0 1
    '''

    atom_name = [
      'O1', 'C6', 'C5', 'C4', 'C3', 'C2', 'H41', 'H31', 'H21', 'H22', 'H51', 'H61'
    ]

    return textwrap.dedent(zmatrix), atom_name

  def monomer_b_aromatic_zmatrix(self):

    zmatrix = '''\
        X21    :1  DISTANCE  :2   ANGLE     :3   90.0000
        O21   X21   1.4950   :1   90.0000   :2  180.0000
        C21   O21   1.3662   X21   59.0000  :1   90.0000
        C22   C21   1.3665   O21  122.6854   X21   DIHEDRAL
        C23   C22   1.4610   C21  117.7425   O21    8.3001
        C24   C23   1.3618   C22  118.3221   C21  -19.0238
        C25   C24   1.5059   C23  116.1256   C22   -7.3457
        H21   C23   1.0938   C22  119.8205   C21  155.4737
        H22   C24   1.0931   C23  122.9298   C22  173.7068
        H23   C25   1.0974   O21  105.0592   C21 -174.6970
        H24   C25   1.1101   O21  106.9039   C21   68.8999
        H25   C22   1.0910   C21  119.1790   O21  178.6725
        H26   C21   1.0927   C22  124.3425   C23 -165.4041
        0 1
    '''

    atom_name = [
      'O1', 'C6', 'C5', 'C4', 'C3', 'C2', 'H41', 'H31', 'H21', 'H22', 'H51', 'H61'
    ]

    return textwrap.dedent(zmatrix), atom_name
