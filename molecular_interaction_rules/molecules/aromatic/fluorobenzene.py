#!/usr/bin/env python3
#
# Lennard-Jones-Drill-2: FluoroBenzene
# ------------------------------------

# Imports
# -------
import textwrap

class FluoroBenzene(object):

    '''

    RESI FLUB           0.00     ! C6H5F, fluorobenzene, adm jr.
    GROUP
    ATOM C1   CG2R61  -0.10
    ATOM H1   HGR62    0.15  !      H5     H4
    ATOM C2   CG2R61  -0.115 !       \ ___ /
    ATOM H2   HGR61    0.115 !       C5---C4
    ATOM C3   CG2R61  -0.115 !       /     \
    ATOM H3   HGR61    0.115 ! F6--C6      C3--H3
    ATOM C4   CG2R61  -0.115 !      \\     //
    ATOM H4   HGR61    0.115 !       C1---C2
    ATOM C5   CG2R61  -0.10  !       /     \
    ATOM H5   HGR62    0.15  !      H1     H2
    ATOM C6   CG2R66   0.11
    ATOM F6   FGR1    -0.21

    RESI  fluoro      0.0000
    GROUP
    ATOM C1     CQ2R6A   -0.120   ALPHA   -1.6955  THOLE   0.9629   ! Penalty =  0.013
    ATOM H1     HQR6A     0.123                                     ! Penalty =  0.014
    ATOM C2     CQ2R6A   -0.147   ALPHA   -1.7319  THOLE   0.9558   ! Penalty =  0.010
    ATOM H2     HQR6A     0.125                                     ! Penalty =  0.010
    ATOM C3     CQ2R6A   -0.097   ALPHA   -1.8986  THOLE   0.9887   ! Penalty =  0.007
    ATOM H3     HQR6A     0.105                                     ! Penalty =  0.004
    ATOM C4     CQ2R6A   -0.147   ALPHA   -1.7319  THOLE   0.9558   ! Penalty =  0.010
    ATOM H4     HQR6A     0.125                                     ! Penalty =  0.010
    ATOM C5     CQ2R6A   -0.119   ALPHA   -1.6955  THOLE   0.9629   ! Penalty =  0.013
    ATOM H5     HQR6A     0.123                                     ! Penalty =  0.014
    ATOM C6     CQ2R6A    0.247   ALPHA   -1.4107  THOLE   0.9392   ! Penalty =  0.024
    ATOM F      FQR1     -0.218   ALPHA   -0.7689  THOLE   1.1151   ! Penalty =  0.074

    Rule 1,2

    '''

    __CGENFF_ATOM_TYPES__ = {
      'RC1': ['CG2R61', 'CG2R66'],
      'F1': ['FGR1'],
      'H1': ['HGR62', 'CG2R61']
    }

    __DGENFF_ATOM_TYPES__ = {
    }

    def __init__(self):

        self.resi_name = 'FLUB'

    def get_monomer_a_species(self):

        '''

        Get the Monomer A Species

        '''

        monomer_a_species = {
            'RC1': self.monomer_a_aromatic_zmatrix(),
            'F1': self.monomer_a_halogen_zmatrix(),
            'H1': self.get_monomer_a_ortho_hydrogen()
        }

        return monomer_a_species

    def get_monomer_b_species(self):

        monomer_b_species = {
            'F1': self.monomer_b_halogen_zmatrix(),
            'RC1': self.monomer_b_pi_stack_zmatrix(),
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
            F11  C11  1.3688 C12  120.0000 C13  180.0000
            H11  C12  1.0756 C11  120.0000 C13  180.0000
            H12  C13  1.0756 C12  120.0000 C11  180.0000
            H13  C14  1.0756 C13  120.0000 C12  180.0000
            H14  C15  1.0756 C14  120.0000 C13  180.0000
            H15  C16  1.0756 C15  120.0000 C11  180.0000
            0 1
        '''

        atom_name = [
            'C6', 'C5', 'C4', 'C3', 'C2', 'C1',  'F6', 'H5', 'H4', 'H3', 'H2', 'H1'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def monomer_a_halogen_zmatrix(self):

        zmatrix = '''\
            F11
            C11  F11  1.3688
            C12  C11  1.3986   F11  120.0000
            C13  C12  1.3774   C11  120.0000   F11  -180.0000
            C14  C13  1.3774   C12  120.0000   C11     0.0000
            C15  C14  1.3774   C13  120.0000   C12     0.0000
            C16  C11  1.3774   C12  120.0000   C13    -0.0000
            H11  C13  1.0756   C12  120.0000   C11   180.0000
            H12  C14  1.0756   C13  120.0000   C12  -180.0000
            H13  C15  1.0756   C14  120.0000   C13   180.0000
            H14  C16  1.0756   C11  120.0000   F11    -0.0000
            H15  C12  1.0756   C11  120.0000   F11     0.0000
            X11  F11  1.0000   C11   90.0000   C12   180.0000
            0 1
        '''
        atom_name = [
            'F6', 'C6', 'C5', 'C4', 'C3', 'C2', 'C1', 'H5', 'H4', 'H3', 'H2', 'H1'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_a_ortho_hydrogen(self):

        zmatrix = '''\
          H11
          C12  H11  1.0756
          C11  C12  1.3986   H11  120.0000
          F11  C11  1.3688   C12  120.0000   H11     0.0000
          C13  C12  1.3774   C11  120.0000   F11  -180.0000
          C14  C13  1.3774   C12  120.0000   C11     0.0000
          C15  C14  1.3774   C13  120.0000   C12     0.0000
          C16  C11  1.3774   C12  120.0000   C13    -0.0000
          H12  C13  1.0756   C12  120.0000   C11   180.0000
          H13  C14  1.0756   C13  120.0000   C12  -180.0000
          H14  C15  1.0756   C14  120.0000   C13   180.0000
          H15  C16  1.0756   C11  120.0000   F11    -0.0000
          X11  H11  1.0000   C12   90.0000   C11   180.0000
          0 1
        '''

        atom_name = [
          'H1', 'C5', 'C6', 'F6', 'C4', 'C3', 'C2', 'C1', 'H5', 'H4', 'H3', 'H2',
        ]

        return textwrap.dedent(zmatrix), atom_name

    def monomer_b_halogen_zmatrix(self):

        zmatrix = '''\
            F21  :1  DISTANCE   :2   ANGLE     :3   90.0000
            X21  F21  1.0000    :1  90.0000    :2    0.0000
            C21  F21  1.3688    X21  90.0000   :1  DIHEDRAL
            C22  C21  1.3986   F21  120.0000   :1    0.0000
            C23  C22  1.3774   C21  120.0000   F21  -180.0000
            C24  C23  1.3774   C22  120.0000   C21     0.0000
            C25  C24  1.3774   C23  120.0000   C22     0.0000
            C26  C21  1.3774   C22  120.0000   C23    -0.0000
            H21  C23  1.0756   C22  120.0000   C21   180.0000
            H22  C24  1.0756   C23  120.0000   C22  -180.0000
            H23  C25  1.0756   C24  120.0000   C23   180.0000
            H24  C26  1.0756   C21  120.0000   F21    -0.0000
            H25  C22  1.0756   C21  120.0000   F21     0.0000
            0 1
        '''

        atom_name = [
           'F6', 'C6', 'C5', 'C4', 'C3', 'C2', 'C1', 'H5', 'H4', 'H3', 'H2', 'H1'
        ],

        return textwrap.dedent(zmatrix), atom_name

    def monomer_b_pi_stack_zmatrix(self):

        zmatrix = '''\
            X21   :1  DISTANCE  :2   ANGLE     :3   90.0000
            C21  X21  1.3940    :1   90.0000   :2  0.0000
            C22  C21  1.3774 X21   60.0000     :1   DIHEDRAL
            C23  C22  1.3774 C21  120.0000 X21    0.0000
            C24  C23  1.3774 C22  120.0000 C21    0.0000
            C25  C24  1.3774 C23  120.0000 C22    0.0000
            C26  C25  1.3774 C24  120.0000 C23    0.0000
            F21  C21  1.3688 C22  120.0000 C23  180.0000
            H21  C22  1.0756 C21  120.0000 C23  180.0000
            H22  C23  1.0756 C22  120.0000 C21  180.0000
            H23  C24  1.0756 C23  120.0000 C22  180.0000
            H24  C25  1.0756 C24  120.0000 C23  180.0000
            H25  C26  1.0756 C25  120.0000 C21  180.0000
            0 1
          '''

        atom_name = [
          'C6', 'C5', 'C4', 'C3', 'C2', 'C1',  'F6', 'H5', 'H4', 'H3', 'H2', 'H1'
        ],

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_b_ortho_hydrogen(self):

        zmatrix = '''\
            H21    :1  DISTANCE  :2   ANGLE    :3    90.00000
            X21   H21  1.0000    :1   90.0000    :2    0.0000
            C22  H21  1.0756   X21   90.0000     :1  DIHEDRAL
            C21  C22  1.3986   H21  120.0000     :1    0.0000
            F21  C21  1.3688   C22  120.0000   H21     0.0000
            C23  C22  1.3774   C21  120.0000   F21  -180.0000
            C24  C23  1.3774   C22  120.0000   C21     0.0000
            C25  C24  1.3774   C23  120.0000   C22     0.0000
            C26  C21  1.3774   C22  120.0000   C23    -0.0000
            H22  C23  1.0756   C22  120.0000   C21   180.0000
            H23  C24  1.0756   C23  120.0000   C22  -180.0000
            H24  C25  1.0756   C24  120.0000   C23   180.0000
            H25  C26  1.0756   C21  120.0000   F21    -0.0000
            0 1
          '''

        atom_name = [
          'H1', 'C5', 'C6', 'F6', 'C4', 'C3', 'C2', 'C1', 'H5', 'H4', 'H3', 'H2',
        ]

        return textwrap.dedent(zmatrix), atom_name
