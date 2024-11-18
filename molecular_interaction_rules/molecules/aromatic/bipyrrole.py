#!/usr/bin/env python3
#
# Lennard-Jones-Drill-2: Bipyrrole
# --------------------------------

# Imports
# -------
import textwrap


class Bipyrrole(object):

    '''

    RESI 13BPO        0.00 ! C8H8N2, 1,3-bipyrrole, lf

    GROUP
    ATOM N1     NG2R57 -0.35
    ATOM C2     CG2R51 -0.04
    ATOM H2     HGR52   0.14
    ATOM C3     CG2R51 -0.25
    ATOM H3     HGR51   0.15
    ATOM C4     CG2R51 -0.25 ! H4     H5    H9      H10
    ATOM H4     HGR51   0.15 !  \     /      \      /
    ATOM C5     CG2R51 -0.04 !   C4==C5       C9==C10
    ATOM H5     HGR52   0.14 !   |     \     /     |
    ATOM N6     NG2R51 -0.35 !   |      N1--C8     |
    ATOM H6     HGP1    0.35 !   |     /     \     |
    ATOM C7     CG2R51 -0.04 !   C3==C2       C7--N6
    ATOM H7     HGR52   0.14 !  /     \      /      \
    ATOM C8     CG2R57  0.25 ! H3     H2    H7      H6
    ATOM C9     CG2R51 -0.25
    ATOM H9     HGR51   0.15
    ATOM C10    CG2R51 -0.04
    ATOM H10    HGR52   0.14

    RESI  ACSF      0.0000
    GROUP
    ATOM C1     CQ2R5A   -0.114   ALPHA   -1.5686  THOLE   1.0929   ! Penalty =  1.800
    ATOM C2     CQ2R5A   -0.217   ALPHA   -1.7779  THOLE   1.1203   ! Penalty =  0.708
    ATOM C3     CQ2R5A   -0.217   ALPHA   -1.7779  THOLE   1.1203   ! Penalty =  0.708
    ATOM C4     CQ2R5A   -0.115   ALPHA   -1.5686  THOLE   1.0929   ! Penalty =  1.800
    ATOM N1     NQ2R5H    0.103   ALPHA   -1.6699  THOLE   1.0549   ! Penalty = 10.896
    ATOM C5     CQ2R5G    0.041   ALPHA   -1.7648  THOLE   1.2355   ! Penalty =  2.144
    ATOM C6     CQ2R5A   -0.157   ALPHA   -1.5330  THOLE   1.1709   ! Penalty =  0.591
    ATOM N2     NQ2R5A   -0.077   ALPHA   -1.6228  THOLE   1.0473   ! Penalty =  3.201
    ATOM C7     CQ2R5A   -0.222   ALPHA   -1.7127  THOLE   1.0790   ! Penalty =  0.439
    ATOM C8     CQ2R5A   -0.148   ALPHA   -1.6922  THOLE   1.1562   ! Penalty =  0.591
    ATOM H1     HQP1A     0.227                                     ! Penalty =  0.932
    ATOM H2     HQR5A     0.119                                     ! Penalty =  2.386
    ATOM H3     HQP1C     0.128                                     ! Penalty =  0.613
    ATOM H4     HQP1C     0.128                                     ! Penalty =  0.613
    ATOM H5     HQR5A     0.119                                     ! Penalty =  2.386
    ATOM H6     HQR5A     0.127                                     ! Penalty =  5.920
    ATOM H7     HQR5A     0.163                                     ! Penalty =  0.933
    ATOM H8     HQP1C     0.112                                     ! Penalty =  0.563

    Rule 1,2

    '''

    # Covered

    __CGENFF_ATOM_TYPES__ = {
      'RC1': ['NG2R51', 'CG2R51', 'CG2R57'],
      'H1': ['HGP1', 'NG2R57']
    }

    __DGENFF_ATOM_TYPES__ = {
      'RC1': ['NQ2R5A', 'CQ2R5A', 'CQ2R5G'],
      'H1': ['HQP1A', 'NQ2R5H']
    }

    def __init__(self):

        self.resi_name = '13BPO'

    def get_monomer_a_species(self):

        '''

        Get the Monomer A Species

        '''

        monomer_a_species = {
          'RC1': self.get_monomer_a_aromatic_zmatrix(),
          'H1': self.get_monomer_a_hydrogen_zmatrix()
        }

        return monomer_a_species

    def get_monomer_b_species(self):

        monomer_b_species = {
            'RC1': self.get_monomer_b_aromatic_zmatrix(),
            'N1': self.get_monomer_b_hydrogen_zmatrix()
        }

        return monomer_b_species

    def get_monomer_a_aromatic_zmatrix(self):

        zmatrix = '''\
          X11
          N11   X11  1.1000
          C11   N11  1.3766  X11   59.0000
          C12   C11  1.3969  N11  107.6333  X11    0.0000
          C13   C12  1.4274  C11  106.9357  N11   -0.0000
          C14   N11  1.3794  C11  110.7158  C12   -0.0000
          N12   C13  1.4129  C14  125.3974  N11 -180.0000
          C15   N12  1.3834  C13  125.3974  C12   -0.0000
          C16   C15  1.3973  N12  108.0174  C13 -180.0000
          C17   C15  2.2753  N12   71.1992  C13 -180.0000
          C18   N12  1.3831  C13  125.3974  C12 -180.0000
          H11   C14  1.0863  N11  121.5865  C11 -180.0000
          H12   C15  1.0863  N12  121.5865  C13    0.0000
          H13   C16  1.0863  C15  125.3974  N12 -180.0000
          H14   C17  1.0863  C16  127.2862  C15 -180.0000
          H15   C18  1.0863  N12  121.5865  C13   -0.0000
          H16   C12  1.0863  C11  125.3974  N11 -180.0000
          H17   C11  1.0863  C12  130.8095  C13 -180.0000
          H18   N11  1.0120  C11  124.9345  C12 -180.0000
          0 1
        '''

        atom_name = [
          'N6', 'C10', 'C9', 'C8', 'C7', 'N1', 'C5', 'C4', 'C3', 'C2', 'H7', 'H8', 'H9', 'H3', 'H2', 'H10', 'H6'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_a_hydrogen_zmatrix(self):

      zmatrix = '''\
          H18
          N11   H18  1.0120
          C11   N11  1.3766  H18  124.9345
          C12   C11  1.3969  N11  107.6333  H18  -180.0000
          C13   C12  1.4274  C11  106.9357  N11   -0.0000
          C14   N11  1.3794  C11  110.7158  C12   -0.0000
          N12   C13  1.4129  C14  125.3974  N11 -180.0000
          C15   N12  1.3834  C13  125.3974  C12   -0.0000
          C16   C15  1.3973  N12  108.0174  C13 -180.0000
          C17   C15  2.2753  N12   71.1992  C13 -180.0000
          C18   N12  1.3831  C13  125.3974  C12 -180.0000
          H11   C14  1.0863  N11  121.5865  C11 -180.0000
          H12   C15  1.0863  N12  121.5865  C13    0.0000
          H13   C16  1.0863  C15  125.3974  N12 -180.0000
          H14   C17  1.0863  C16  127.2862  C15 -180.0000
          H15   C18  1.0863  N12  121.5865  C13   -0.0000
          H16   C12  1.0863  C11  125.3974  N11 -180.0000
          H17   C11  1.0863  C12  130.8095  C13 -180.0000
          X11   H18  1.0000  N11   90.0000  C11    0.0000
          0 1
        '''

      atom_name = [
        'H6', 'N6', 'C10', 'C9', 'C8', 'C7', 'N1', 'C5', 'C4', 'C3', 'C2', 'H7', 'H8', 'H9', 'H3', 'H2', 'H10',
      ]

      return textwrap.dedent(zmatrix), atom_name

    def get_monomer_b_aromatic_zmatrix(self):

      zmatrix = '''\
          X21    :1  DISTANCE  :2   ANGLE   :3    90.0000
          N21   X21  1.1000    :1   90.0000    :2  180.000
          C21   N21  1.3766  X21   59.0000     :1   90.0000
          C22   C21  1.3969  N21  107.6333  X21   DIHEDRAL
          C23   C22  1.4274  C21  106.9357  N21   -0.0000
          C24   N21  1.3794  C21  110.7158  C22   -0.0000
          N22   C23  1.4129  C24  125.3974  N21 -180.0000
          C25   N22  1.3834  C23  125.3974  C22   -0.0000
          C26   C25  1.3973  N22  108.0174  C23 -180.0000
          C27   C25  2.2753  N22   71.1992  C23 -180.0000
          C28   N22  1.3831  C23  125.3974  C22 -180.0000
          H21   C24  1.0863  N21  121.5865  C21 -180.0000
          H22   C25  1.0863  N22  121.5865  C23    0.0000
          H23   C26  1.0863  C25  125.3974  N22 -180.0000
          H24   C27  1.0863  C26  127.2862  C25 -180.0000
          H25   C28  1.0863  N22  121.5865  C23   -0.0000
          H26   C22  1.0863  C21  125.3974  N21 -180.0000
          H27   C21  1.0863  C22  130.8095  C23 -180.0000
          H28   N21  1.0120  C21  124.9345  C22 -180.0000
          0 1
        '''

      atom_name = [
        'N6', 'C10', 'C9', 'C8', 'C7', 'N1', 'C5', 'C4', 'C3', 'C2', 'H7', 'H8', 'H9', 'H3', 'H2', 'H10', 'H6'
      ]

      return textwrap.dedent(zmatrix), atom_name

    def get_monomer_b_hydrogen_zmatrix(self):

      zmatrix = '''\
          H28   :1  DISTANCE    :2   ANGLE    :3   90.0000
          X21   H28  1.0000     :1   90.0000  :2    0.0000
          N21   H28  1.0120     X21   90.0000 :1  DIHEDRAL
          C21   N21  1.3766  H28  124.9345    :1    0.0000
          C22   C21  1.3969  N21  107.6333  H28  -180.0000
          C23   C22  1.4274  C21  106.9357  N21   -0.0000
          C24   N21  1.3794  C21  110.7158  C22   -0.0000
          N22   C23  1.4129  C24  125.3974  N21 -180.0000
          C25   N22  1.3834  C23  125.3974  C22   -0.0000
          C26   C25  1.3973  N22  108.0174  C23 -180.0000
          C27   C25  2.2753  N22   71.1992  C23 -180.0000
          C28   N22  1.3831  C23  125.3974  C22 -180.0000
          H21   C24  1.0863  N21  121.5865  C21 -180.0000
          H22   C25  1.0863  N22  121.5865  C23    0.0000
          H23   C26  1.0863  C25  125.3974  N22 -180.0000
          H24   C27  1.0863  C26  127.2862  C25 -180.0000
          H25   C28  1.0863  N22  121.5865  C23   -0.0000
          H26   C22  1.0863  C21  125.3974  N21 -180.0000
          H27   C21  1.0863  C22  130.8095  C23 -180.0000
          0 1
        '''

      atom_name = [
        'H6', 'N6', 'C10', 'C9', 'C8', 'C7', 'N1', 'C5', 'C4', 'C3', 'C2', 'H7', 'H8', 'H9', 'H3', 'H2', 'H10'
      ]

      return textwrap.dedent(zmatrix), atom_name


