#!/usr/bin/env python3
#
# Lennard-Jones-Drill-2: DimethylSulfone
# --------------------------------------

# Imports
# -------

import textwrap

class DimethylSulfone(object):

    '''

    RESI DMSN          0.00 ! C2H6O2S, dimethyl sulfone, xhe
    GROUP
    ATOM S    SG3O2    0.14
    ATOM O1   OG2P1   -0.36 !     H31
    ATOM O2   OG2P1   -0.36 !      |
    ATOM C3   CG331    0.02 !  H32-C3-H33
    ATOM H31  HGA3     0.09 !      |
    ATOM H32  HGA3     0.09 !   O1=S=O2
    ATOM H33  HGA3     0.09 !      |
    ATOM C4   CG331    0.02 !  H42-C4-H43
    ATOM H41  HGA3     0.09 !      |
    ATOM H42  HGA3     0.09 !     H41
    ATOM H43  HGA3     0.09

    RESI  methyl      0.0000
    GROUP
    ATOM C1     CQ33A    -0.281   ALPHA   -1.8726  THOLE   0.9529   ! Penalty =  0.430
    ATOM H1     HQA3A     0.109                                     ! Penalty =  0.105
    ATOM H2     HQA3A     0.109                                     ! Penalty =  0.105
    ATOM H3     HQA3A     0.109                                     ! Penalty =  0.105
    ATOM S      SQ3OB     0.985   ALPHA   -2.0445  THOLE   0.9153   ! Penalty =  0.530
    ATOM C2     CQ33A    -0.280   ALPHA   -1.8726  THOLE   0.9529   ! Penalty =  0.430
    ATOM H4     HQA3A     0.109                                     ! Penalty =  0.105
    ATOM H5     HQA3A     0.109                                     ! Penalty =  0.105
    ATOM H6     HQA3A     0.109                                     ! Penalty =  0.105
    ATOM O1     OQ2C2B   -0.539   ALPHA   -1.1555  THOLE   0.9666   ! Penalty =  0.199
    ATOM O2     OQ2C2B   -0.539   ALPHA   -1.1555  THOLE   0.9666   ! Penalty =  0.199

    Rule 3
    '''

    __CGENFF_ATOM_TYPES__ = {
        'S1': ['SG3O2'],
        'O1': ['OG2P1']
    }

    __DGENFF_ATOM_TYPES__ = {
      'S1': ['SQ3OB'],
      'O1': ['OQ2C2B']
    }

    def __init__(self):

        self.resi_name = 'DMSN'

    def get_monomer_a_species(self):

        '''

        Get the Monomer A Species

        '''

        monomer_a_species = {
            'O1': self.get_oxygen_zmatrix(),
            'S1': self.get_sulphur_matrix()
        }

        return monomer_a_species

    def get_oxygen_zmatrix(self):

        zmatrix = '''\
          O11
          S11 O11 1.4948
          C11 S11 1.7995 O11  107.8073
          H11 C11 1.0984 S11  108.6706 O11 -175.3736
          H12 C11 1.0984 S11  108.6706 O11  -52.7697
          H13 C11 1.0984 S11  105.2498 O11   65.9283
          C12 S11 1.7995 C11  103.5586 H11  -61.3019
          H14 C12 1.0984 S11  108.6706 O11   52.7696
          H15 C12 1.0984 S11  105.2498 O11  -65.9283
          H16 C12 1.0984 S11  108.6706 O11  175.3736
          O12 S11 1.4948 C11  107.8074 H11   52.7696
          X11 O11 1.0000 S11   90.0000 C11  180.0000
          0 1
        '''

        atom_name = [
          'O1', 'S', 'C3', 'H31', 'H32', 'H33', 'C4', 'H41', 'H42', 'H43', 'O2'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_sulphur_matrix(self):

        zmatrix = '''\
          S11
          O11 S11 1.4948
          C11 S11 1.7995 O11  107.8073
          H11 C11 1.0984 S11  108.6706 O11 -175.3736
          H12 C11 1.0984 S11  108.6706 O11  -52.7697
          H13 C11 1.0984 S11  105.2498 O11   65.9283
          C12 S11 1.7995 C11  103.5586 H11  -61.3019
          H14 C12 1.0984 S11  108.6706 O11   52.7696
          H15 C12 1.0984 S11  105.2498 O11  -65.9283
          H16 C12 1.0984 S11  108.6706 O11  175.3736
          O12 S11 1.4948 C11  107.8074 H11   52.7696
          X11 S11 1.0000 O11  55.0000 C11  135.0000
          0 1
        '''

        atom_name = [
          'S', 'O1','C3', 'H31', 'H32', 'H33', 'C4', 'H41', 'H42', 'H43', 'O2'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_b_species(self):

        '''

        Get the Monomer B Species

        '''

        monomer_b_species = {
            'O1': self.get_monomer_b_oxygen_zmatrix(),
            'S1': self.get_monomer_b_sulphur_zmatrix()
        }

        return monomer_b_species

    def get_monomer_b_oxygen_zmatrix(self):

        zmatrix = '''\
            O21  :1  DISTANCE   :2  ANGLE    :3    DIHEDRAL
            S21 O21 1.4948      :1  180.0000  :2   180.0000
            C21 S21 1.7995 O21  107.8073      :1     0.0000
            H21 C21 1.0984 S21  108.6706 O21 -175.3736
            H22 C21 1.0984 S21  108.6706 O21  -52.7697
            H23 C21 1.0984 S21  105.2498 O21   65.9283
            C22 S21 1.7995 C21  103.5586 H21  -61.3019
            H24 C22 1.0984 S21  108.6706 O21   52.7696
            H25 C22 1.0984 S21  105.2498 O21  -65.9283
            H26 C22 1.0984 S21  108.6706 O21  175.3736
            O22 S21 1.4948 C21  107.8074 H21   52.7696
            X21 O21 1.0000 S21   90.0000 C21  180.0000
            0 1
          '''

        atom_name = [
          'O1', 'S', 'C3', 'H31', 'H32', 'H33', 'C4', 'H41', 'H42', 'H43', 'O2'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_b_sulphur_zmatrix(self):

        zmatrix = '''\
              S21  :1  DISTANCE   :2  ANGLE    :3    DIHEDRAL
              O21 S21 1.4948      :1  180.0000  :2   180.0000
              C21 S21 1.7995 O21  107.8073      :1     0.0000
              H21 C21 1.0984 S21  108.6706 O21 -175.3736
              H22 C21 1.0984 S21  108.6706 O21  -52.7697
              H23 C21 1.0984 S21  105.2498 O21   65.9283
              C22 S21 1.7995 C21  103.5586 H21  -61.3019
              H24 C22 1.0984 S21  108.6706 O21   52.7696
              H25 C22 1.0984 S21  105.2498 O21  -65.9283
              H26 C22 1.0984 S21  108.6706 O21  175.3736
              O22 S21 1.4948 C21  107.8074 H21   52.7696
              X21 O21 1.0000 S21   90.0000 C21  180.0000
              0 1
            '''

        atom_name = [
          'S', 'O1', 'C3', 'H31', 'H32', 'H33', 'C4', 'H41', 'H42', 'H43', 'O2'
        ]

        return textwrap.dedent(zmatrix), atom_name
