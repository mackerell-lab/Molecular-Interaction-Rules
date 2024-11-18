#!/usr/bin/env python3
#
# Lennard-Jones-Drill-2: Methyleneoxindole
# ----------------------------------------

# Imports
# -------
import textwrap

class Methyleneoxindole(object):

  '''

  RESI MEOI         0.000 ! C9H7NO, methyleneoxindole, kevo & xxwy

  GROUP
  ATOM C5   CG2DC3  -0.41
  ATOM H51  HGA5     0.21
  ATOM H52  HGA5     0.21
  ATOM C6   CG25C1  -0.10 !            H51    H52
  ATOM C7   CG2R53   0.67 !               \  /
  ATOM O7   OG2D1   -0.57 !         H13    C5
  ATOM N8   NG2R51  -0.45 !          |     ||
  ATOM H8   HGP1     0.32 !         C13    C6   O7
  ATOM C9   CG2RC0   0.25 !        //  \  /  \ //
  ATOM C10  CG2R61  -0.34 ! H12--C12    C14   C7
  ATOM H10  HGR61    0.24 !       |     ||    |
  ATOM C11  CG2R61  -0.20 ! H11--C11    C9----N8
  ATOM H11  HGR61    0.22 !        \\  / \
  ATOM C12  CG2R61  -0.23 !         C10  H8
  ATOM H12  HGR61    0.21 !          |
  ATOM C13  CG2R61  -0.32 !         H10
  ATOM H13  HGR61    0.28
  ATOM C14  CG2RC0   0.01

  RESI  3methy      0.0000
  GROUP
  ATOM O      OQ2C1A    0.001   ALPHA   -1.2596  THOLE   0.6736   ! Penalty =  1.475
  ATOM N      NQ2R5A   -0.217   ALPHA   -1.5134  THOLE   1.2194   ! Penalty = 12.037
  ATOM C1     CQ2R6F    0.061   ALPHA   -1.6095  THOLE   1.5389   ! Penalty = 100.679
  ATOM C2     CQ2R6F    0.025   ALPHA   -1.5217  THOLE   1.7787   ! Penalty =  1.449
  ATOM C3     CQ25C1    0.029   ALPHA   -1.7104  THOLE   1.2187   ! Penalty =  7.317
  ATOM C4     CQ2R5B    0.362   ALPHA   -1.3774  THOLE   1.2479   ! Penalty =  3.413
  ATOM C5     CQ2R6A   -0.122   ALPHA   -1.7625  THOLE   1.2030   ! Penalty =  0.498
  ATOM C6     CQ2R6A   -0.241   ALPHA   -1.7300  THOLE   1.2387   ! Penalty =  0.233
  ATOM C7     CQ2R6A   -0.122   ALPHA   -1.8197  THOLE   1.1777   ! Penalty =  0.172
  ATOM C8     CQ2R6A   -0.153   ALPHA   -1.8521  THOLE   1.2886   ! Penalty =  0.065
  ATOM C9     CQ2DC3   -0.383   ALPHA   -1.6137  THOLE   1.1132   ! Penalty = 75.325
  ATOM H1     HQP1A     0.327                                     ! Penalty =  0.516
  ATOM H2     HQR6A     0.119                                     ! Penalty =  0.396
  ATOM H3     HQR6A     0.142                                     ! Penalty =  0.132
  ATOM H4     HQR6A     0.116                                     ! Penalty =  0.088
  ATOM H5     HQR6A     0.121                                     ! Penalty =  0.017
  ATOM H6     HQ2C1B    0.164                                     ! Penalty = 41.892
  ATOM H7     HQ2C1B    0.164                                     ! Penalty = 41.892
  ATOM LP1    LPQO1    -0.228                                     ! Penalty =  0.186
  ATOM LP2    LPQO1    -0.165                                     ! Penalty =  0.186

  Rule 2

  '''

  # Covered

  __CGENFF_ATOM_TYPES__ = {
    'RC1': ['CG25C1', 'CG2RC0'],
  }

  __DGENFF_ATOM_TYPES__ = {
    'RC1': ['CQ2R6A', 'CQ2R6F'],
  }

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



