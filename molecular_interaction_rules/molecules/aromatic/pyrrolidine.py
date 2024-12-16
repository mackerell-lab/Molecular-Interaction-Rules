#!/usr/bin/env python3
#
# Lennard-Jones-Drill-2: Pyrrolidine
# ----------------------------------

# Imports
# -------
import textwrap

class Pyrrolidine(object):

  def __init__(self):

      self.resi_name = 'PRLD'

  def get_monomer_a_species(self):

      monomer_a_species = {
        'H1': self.get_monomer_a_hydrogen_zmatrix(),
        'RC1': self.get_monomer_a_aromatic_zmatrix()
      }

      return monomer_a_species

  def get_monomer_b_species(self):

      monomer_b_species = {
        'RC1': self.get_monomer_b_aromatic_zmatrix(),
        'H1': self.get_monomer_b_hydrogen_zmatrix()
      }

      return monomer_b_species

  def get_monomer_a_hydrogen_zmatrix(self):

      zmatrix = '''\
        H11
        N11   H11    1.0243
        C11   N11    1.4776    H11  107.9129
        C12   C11    1.5529    N11  106.8472   H11  -75.7702
        C13   C12    1.5579    C11  104.0259   N11  -23.3895
        C14   N11    1.4774    C11  102.8883   C12   38.1662
        H12   C13    1.1000    C12  112.4064   C11 -121.2645
        H13   C13    1.1000    C12  110.5580   C11  118.4384
        H14   C14    1.1000    N11  107.4234   C11   80.0279
        H15   C14    1.1000    N11  110.5580   C11 -162.3234
        H16   C12    1.1000    C11  110.5580   N11   95.3232
        H17   C12    1.1000    C11  111.9077   N11 -144.9858
        H18   C11    1.1000    C12  113.6638   C13 -145.6173
        H19   C11    1.1000    C12  110.5580   C13   92.9993
        X11   H11    1.0000    N11   90.0000   C11  180.0000
        0 1
      '''

      atom_name = [
        'H1', 'N1', 'C2', 'C3', 'C4', 'C5', 'H42', 'H41', 'H51', 'H52', 'H32', 'H31', 'H22', 'H21'
      ]

      return textwrap.dedent(zmatrix), atom_name

  def get_monomer_a_aromatic_zmatrix(self):

      zmatrix = '''\
        X11
        N11  X11  1.1000
        C11  N11  1.4776  X11   60.0000
        C12  C11  1.5529  N11  106.8472  X11    0.0000
        C13  C12  1.5579  C11  104.0259  N11  -23.3895
        H11  C13  1.1000  C12  112.4064  C11 -121.2645
        H12  C13  1.1000  C12  110.5580  C11  118.4384
        C14  N11  1.4774  C11  102.8883  C12   38.1662
        H12  C14  1.1000  N11  107.4234  C11   80.0279
        H13  C14  1.1000  N11  110.5580  C11 -162.3234
        H14  C12  1.1000  C11  110.5580  N11   95.3232
        H15  C12  1.1000  C11  111.9077  N11 -144.9858
        H16  C11  1.1000  C12  113.6638  C13 -145.6173
        H17  C11  1.1000  C12  110.5580  C13   92.9993
        H18  N11  1.0243  C11  107.9129  C12  -75.7702
        0 1
      '''

      atom_name = [
         'N1', 'C2', 'C3', 'C4', 'C5', 'H42', 'H41', 'H51', 'H52', 'H32', 'H31', 'H22', 'H21','H1',
      ]

      return textwrap.dedent(zmatrix), atom_name

  def get_monomer_b_aromatic_zmatrix(self):

      zmatrix = '''\
         X21   :1  DISTANCE   :2  ANGLE    :3   90.0000
        N21  X21  1.1000     :1  90.0000    :2  180.0000
        C21  N21  1.4776  X21   60.0000   :1   90.0000
        C22  C21  1.5529  N21  106.8472  X21   DIHEDRAL
        C23  C22  1.5579  C21  104.0259  N21  -23.3895
        H21  C23  1.1000  C22  112.4064  C21 -121.2645
        H22  C23  1.1000  C22  110.5580  C21  118.4384
        C24  N21  1.4774  C21  102.8883  C22   38.1662
        H22  C24  1.1000  N21  107.4234  C21   80.0279
        H23  C24  1.1000  N21  110.5580  C21 -162.3234
        H24  C22  1.1000  C21  110.5580  N21   95.3232
        H25  C22  1.1000  C21  111.9077  N21 -144.9858
        H26  C21  1.1000  C22  113.6638  C23 -145.6173
        H27  C21  1.1000  C22  110.5580  C23   92.9993
        H28  N21  1.0243  C21  107.9129  C22  -75.7702
        0 1
      '''

      atom_name = [
        'N1', 'C2', 'C3', 'C4', 'C5', 'H42', 'H41', 'H51', 'H52', 'H32', 'H31', 'H22', 'H21','H1',
      ]

      return textwrap.dedent(zmatrix), atom_name

  def get_monomer_b_hydrogen_zmatrix(self):

    zmatrix = '''\
        H21    :1  DISTANCE  :2   ANGLE    :3   90.0000
        X21  H21    1.0000   :1   90.0000   :2   0.0000
        N21   H21    1.0243    X21  90.0000    :1  DIHEDRAL
        C21   N21    1.4776    H21  107.9129   :1  0.0000
        C22   C21    1.5529    N21  106.8472   H21  -75.7702
        C23   C22    1.5579    C21  104.0259   N21  -23.3895
        C24   N21    1.4774    C21  102.8883   C22   38.1662
        H22   C23    1.1000    C22  112.4064   C21 -121.2645
        H23   C23    1.1000    C22  110.5580   C21  118.4384
        H24   C24    1.1000    N21  107.4234   C21   80.0279
        H25   C24    1.1000    N21  110.5580   C21 -162.3234
        H26   C22    1.1000    C21  110.5580   N21   95.3232
        H27   C22    1.1000    C21  111.9077   N21 -144.9858
        H28   C21    1.1000    C22  113.6638   C23 -145.6173
        H29   C21    1.1000    C22  110.5580   C23   92.9993
        0 1
      '''

    atom_name = [
      'H1', 'N1', 'C2', 'C3', 'C4', 'C5', 'H42', 'H41', 'H51', 'H52', 'H32', 'H31', 'H22', 'H21'
    ]

    return textwrap.dedent(zmatrix), atom_name


