#!/usr/bin/env python3
#
# Lennard-Jones-Drill-2: BromoBenzene
# -----------------------------------

# Imports
# -------
import textwrap

class BromoBenzene(object):

  '''

    CGenFF

    RESI BROB          0.000 ! C6H5Br, Bromobenzene + lone pair , isg/fylin
    GROUP                    !
    ATOM C1   CG2R61  -0.06  !
    ATOM H1   HGR62    0.15  !         H5     H4
    ATOM C2   CG2R61  -0.115 !          \ ___ /
    ATOM H2   HGR61    0.115 !          C5---C4
    ATOM C3   CG2R61  -0.115 !          /     \
    ATOM H3   HGR61    0.115 ! LP-BR--C6      C3--H3
    ATOM C4   CG2R61  -0.115 !          \\   //
    ATOM H4   HGR61    0.115 !          C1---C2
    ATOM C5   CG2R61  -0.06  !          / \
    ATOM H5   HGR62    0.15  !         H1     H2
    ATOM C6   CG2R61  -0.05  !
    ATOM BR   BRGR1   -0.18  !
    ATOM LP   LPH      0.05  !

    DGenFF

    RESI  bromob      0.0000
    GROUP
    ATOM C1     CQ2R6A   -0.027   ALPHA   -1.7561  THOLE   0.9915   ! Penalty =  0.520
    ATOM H1     HQR6A     0.083                                     ! Penalty =  0.343
    ATOM C2     CQ2R6A   -0.106   ALPHA   -1.7287  THOLE   1.0120   ! Penalty =  0.381
    ATOM H2     HQR6A     0.099                                     ! Penalty =  0.254
    ATOM C3     CQ2R6A   -0.097   ALPHA   -1.9153  THOLE   0.9977   ! Penalty =  0.192
    ATOM H3     HQR6A     0.091                                     ! Penalty =  0.004
    ATOM C4     CQ2R6A   -0.106   ALPHA   -1.7287  THOLE   1.0120   ! Penalty =  0.381
    ATOM H4     HQR6A     0.099                                     ! Penalty =  0.254
    ATOM C5     CQ2R6A   -0.028   ALPHA   -1.7561  THOLE   0.9915   ! Penalty =  0.520
    ATOM H5     HQR6A     0.083                                     ! Penalty =  0.343
    ATOM C6     CQ2R6A   -0.045   ALPHA   -2.0526  THOLE   0.9329   ! Penalty =  1.145
    ATOM BR     BRQR1    -0.096   ALPHA   -3.2077  THOLE   1.0836   ! Penalty =  2.798
    ATOM LP1    LPD       0.050                                     ! Penalty =  0.000

    Rule 1,2

  '''

  __CGENFF_ATOM_TYPES__ = {
      'BR1': ['BRGR1'],
      'RC1': ['CG2R61'],
      'H1': ['HGR62']
  }

  __DGENFF_ATOM_TYPES__ = {
  }

  def __init__(self):

    self.resi_name = 'BROB'

  def get_monomer_a_species(self):

    '''

    Get the Monomer A Species

    '''

    monomer_a_species = {
        'BR1': self.monomer_a_halogen_zmatrix(),
        'RC1': self.monomer_a_aromatic_zmatrix(),
        'H1': self.get_monomer_a_ortho_hydrogen()
    }

    return monomer_a_species

  def get_monomer_b_species(self):

    monomer_b_species = {
        'BR1': self.monomer_b_halogen_zmatrix(),
        'RC1': self.monomer_b_aromatic_zmatrix(),
        'H1': self.get_monomer_b_ortho_hydrogen()
    }

    return monomer_b_species

  def monomer_a_aromatic_zmatrix(self):

      zmatrix = '''\
              X11
              C11  X11  1.3940
              C12  C11  1.3774 X11   60.0000
              C13  C12  1.3774 C11  120.0000 X11    0.0000
              C14  C13  1.3774 C12  120.0000 C11    0.0000
              C15  C14  1.3774 C13  120.0000 C12    0.0000
              C16  C15  1.3774 C14  120.0000 C13    0.0000
              BR11 C11  1.937 C12  120.0000 C13  180.0000
              H11  C12  1.0756 C11  120.0000 C13  180.0000
              H12  C13  1.0756 C12  120.0000 C11  180.0000
              H13  C14  1.0756 C13  120.0000 C12  180.0000
              H14  C15  1.0756 C14  120.0000 C13  180.0000
              H15  C16  1.0756 C15  120.0000 C11  180.0000
              0 1
          '''

      atom_name = [
        'C6', 'C5', 'C4', 'C3', 'C2', 'C1', 'BR', 'H5', 'H4', 'H3', 'H2', 'H1'
      ]

      return textwrap.dedent(zmatrix), atom_name

  def monomer_a_halogen_zmatrix(self):

      zmatrix = '''\
          BR11
          C11  BR11  1.8000
          C12  C11  1.3986   BR11  120.0000
          C13  C12  1.3774   C11  120.0000   BR11  -180.0000
          C14  C13  1.3774   C12  120.0000   C11     0.0000
          C15  C14  1.3774   C13  120.0000   C12     0.0000
          C16  C11  1.3774   C12  120.0000   C13    -0.0000
          H11  C13  1.0756   C12  120.0000   C11   180.0000
          H12  C14  1.0756   C13  120.0000   C12  -180.0000
          H13  C15  1.0756   C14  120.0000   C13   180.0000
          H14  C16  1.0756   C11  120.0000   BR11    -0.0000
          H15  C12  1.0756   C11  120.0000   BR11     0.0000
          X11  BR11  1.0000   C11   90.0000   C12   180.0000
          0 1
      '''
      atom_name = [
        'BR', 'C6', 'C5', 'C4', 'C3', 'C2', 'C1', 'H5', 'H4', 'H3', 'H2', 'H1'
      ]

      return textwrap.dedent(zmatrix), atom_name

  def get_monomer_a_ortho_hydrogen(self):

      zmatrix = '''\
            H11
            C12  H11  1.0756
            C11  C12  1.3986   H11  120.0000
            BR11  C11  1.937   C12  120.0000   H11     0.0000
            C13  C12  1.3774   C11  120.0000   BR11  -180.0000
            C14  C13  1.3774   C12  120.0000   C11     0.0000
            C15  C14  1.3774   C13  120.0000   C12     0.0000
            C16  C11  1.3774   C12  120.0000   C13    -0.0000
            H12  C13  1.0756   C12  120.0000   C11   180.0000
            H13  C14  1.0756   C13  120.0000   C12  -180.0000
            H14  C15  1.0756   C14  120.0000   C13   180.0000
            H15  C16  1.0756   C11  120.0000   BR11    -0.0000
            X11  H11  1.0000   C12   90.0000   C11   180.0000
            0 1
          '''

      atom_name = [
        'H1', 'C5', 'C6', 'BR', 'C4', 'C3', 'C2', 'C1', 'H5', 'H4', 'H3', 'H2',
      ]

      return textwrap.dedent(zmatrix), atom_name

  def monomer_b_halogen_zmatrix(self):

      zmatrix = '''\
         BR21  :1  DISTANCE   :2   ANGLE    :3    90.0000
         X21  BR21  1.0000    :1  90.0000   :2     0.0000
         C21  BR21  1.937    X21  90.0000    :1   DIHEDRAL
         C22  C21  1.3986   BR21  120.0000   :1     0.0000
         C23  C22  1.3774   C21  120.0000   BR21  -180.0000
         C24  C23  1.3774   C22  120.0000   C21     0.0000
         C25  C24  1.3774   C23  120.0000   C22     0.0000
         C26  C21  1.3774   C22  120.0000   C23    -0.0000
         H21  C23  1.0756   C22  120.0000   C21   180.0000
         H22  C24  1.0756   C23  120.0000   C22  -180.0000
         H23  C25  1.0756   C24  120.0000   C23   180.0000
         H24  C26  1.0756   C21  120.0000   BR21    -0.0000
         H25  C22  1.0756   C21  120.0000   BR21     0.0000
         0 1
      '''

      atom_name = [
        'BR', 'C6', 'C5', 'C4', 'C3', 'C2', 'C1', 'H5', 'H4', 'H3', 'H2', 'H1'
      ],

      return textwrap.dedent(zmatrix), atom_name

  def monomer_b_aromatic_zmatrix(self):

      zmatrix = '''\
          X21   :1  DISTANCE  :2   ANGLE   :3   90.0000
          C21  X21  1.3940    :1   90.0000    :2  180.0000
          C22  C21  1.3774 X21   60.0000      :1   90.0000
          C23  C22  1.3774 C21  120.0000 X21    DIHEDRAL
          C24  C23  1.3774 C22  120.0000 C21    0.0000
          C25  C24  1.3774 C23  120.0000 C22    0.0000
          C26  C25  1.3774 C24  120.0000 C23    0.0000
          BR21 C21  1.937  C22  120.0000 C23  180.0000
          H21  C22  1.0756 C21  120.0000 C23  180.0000
          H22  C23  1.0756 C22  120.0000 C21  180.0000
          H23  C24  1.0756 C23  120.0000 C22  180.0000
          H24  C25  1.0756 C24  120.0000 C23  180.0000
          H25  C26  1.0756 C25  120.0000 C21  180.0000
          0 1
      '''

      atom_name = [
        'C6', 'C5', 'C4', 'C3', 'C2', 'C1', 'BR', 'H5', 'H4', 'H3', 'H2', 'H1'
      ]

      return textwrap.dedent(zmatrix), atom_name


  def get_monomer_b_ortho_hydrogen(self):

      zmatrix = '''\
              H21    :1  DISTANCE  :2   ANGLE    :3   90.0000
              X21   H21  1.0000    :1   90.0000    :2    0.0000
              C22  H21  1.0756   X21   90.0000     :1  DIHEDRAL
              C21  C22  1.3986   H21  120.0000     :1    0.0000
              BR21  C21  1.937   C22  120.0000   H21     0.0000
              C23  C22  1.3774   C21  120.0000   BR21  -180.0000
              C24  C23  1.3774   C22  120.0000   C21     0.0000
              C25  C24  1.3774   C23  120.0000   C22     0.0000
              C26  C21  1.3774   C22  120.0000   C23    -0.0000
              H22  C23  1.0756   C22  120.0000   C21   180.0000
              H23  C24  1.0756   C23  120.0000   C22  -180.0000
              H24  C25  1.0756   C24  120.0000   C23   180.0000
              H25  C26  1.0756   C21  120.0000   BR21    -0.0000
              0 1
            '''

      atom_name = [
        'H1', 'C5', 'C6', 'BR', 'C4', 'C3', 'C2', 'C1', 'H5', 'H4', 'H3', 'H2',
      ]

      return textwrap.dedent(zmatrix), atom_name

