#!/usr/bin/env python3
#
# Lennard-Jones-Drill-2: Imidazole
# --------------------------------

# Imports
# -------
import textwrap

class Imidazole(object):

    '''

    CGenFF

    GROUP
    ATOM CG   CG2R51  -0.05 !       HD1    HE1
    ATOM HG   HGR52    0.09 !        |     /
    ATOM CD2  CG2R51   0.22 !       ND1--CE1
    ATOM HD2  HGR52    0.10 !      /     ||
    ATOM ND1  NG2R51  -0.36 ! HG-CG      ||
    ATOM HD1  HGP1     0.32 !     \\     ||
    ATOM CE1  CG2R53   0.25 !      CD2--NE2
    ATOM HE1  HGR52    0.13 !       |
    ATOM NE2  NG2R50  -0.70 !      HD2

    DGenFF

    RESI IMID          0.000 ! Imidazole
    !
    !       HD1    HE1
    !        |     /
    !       ND1--CE1
    !      /     ||
    ! HG-CG      ||
    !     \\     ||
    !      CD2--NE2
    !       |      :
    !      HD2
    !

    GROUP
    ATOM CG   CD2R5A  -0.343  ALPHA -1.632  THOLE 0.949
    ATOM HG   HDR5A    0.215
    ATOM CD2  CD2R5A   0.037  ALPHA -1.626  THOLE 0.949
    ATOM HD2  HDR5A    0.091
    ATOM ND1  ND2R5A  -0.053  ALPHA -1.393  THOLE 0.670
    ATOM HD1  HDP1A    0.244
    ATOM CE1  CD2R5B   0.040  ALPHA -1.632  THOLE 1.568
    ATOM HE1  HDR5B    0.107
    ATOM NE2  ND2R5B   0.000  ALPHA -1.044  THOLE 0.788
    ATOM LP1  LPD     -0.338

    Rule 2, 3
    '''

    # Covered

    __CGENFF_ATOM_TYPES__ = {
      'RC1': ['CG2R51', 'NG2R51', 'CG2R53', 'NG2R51'],
      'N1': ['ND2R5B']
    }

    __DGENFF_ATOM_TYPES__ = {
      'RC1': ['CD2R5A', 'CD2R5B', 'ND2R5B', 'ND2R5A'],
      'N1': ['NG2R50']
    }

    def __init__(self):

        self.resi_name = 'IMID'

    def get_monomer_a_species(self):

        '''

        Get the Monomer A Species

        '''

        monomer_a_species = {
            'RC1': self.monomer_a_aromatic_zmatrix(),
            'N1': self.monomer_a_nitrogen_zmatrix()
        }

        return monomer_a_species

    def get_monomer_b_species(self):

        monomer_b_species = {
            'N1': self.monomer_b_nitrogen_zmatrix(),
            'RC1': self.monomer_b_aromatic_zmatrix()
        }

        return monomer_b_species

    def monomer_a_aromatic_zmatrix(self):

        zmatrix = '''\
            X11
            C11  X11  1.1134
            C12  C11  1.3918  X11   59.0000
            N11  C12  1.3834  C11  110.7737   X11   0.0000
            C13  N11  1.3730  C12  105.1514   C11   -0.0000
            N12  C11  1.3815  C12  104.8582   N11    0.0000
            H11  C13  1.0872  N11  126.0493   C12 -180.0000
            H12  N12  1.0130  C11  126.1042   C12  180.0000
            H13  C12  1.0878  N11  121.5378   C13  180.0000
            H14  C11  1.0868  C12  132.6904   N11 -180.0000
            0 1
        '''

        atom_name = [
            'CG', 'CD2', 'NE2', 'CE1', 'ND1', 'HE1', 'HD1', 'HD2', 'HG'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def monomer_a_nitrogen_zmatrix(self):

        zmatrix = '''\
            N11
            C12  N11  1.3834
            C11  C12  1.3918  N11  110.7737
            C13  N11  1.3351  C12  105.1514   C11   -0.0000
            N12  C11  1.3730  C12  104.8582   N11    0.0000
            H11  C13  1.0872  N11  126.0493   C12 -180.0000
            H12  N12  1.0130  C11  126.1042   C12  180.0000
            H13  C12  1.0878  N11  121.5378   C13  180.0000
            H14  C11  1.0868  C12  132.6904   N11 -180.0000
            X11  N11  1.0000  C12   90.0000   C11  180.0000
            0 1
        '''

        atom_name = [
            'NE2', 'CD2', 'CG',  'CE1', 'ND1', 'HE1', 'HD1', 'HD2', 'HG'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def monomer_b_nitrogen_zmatrix(self):

        zmatrix = '''\
            N21   :1  DISTANCE  :2  ANGLE      :3   DIHEDRAL
            C22  N21  1.3730    :1  126.5480   :2     0.0000
            C23  N21  1.3730    :1  126.9043   :2   180.0000
            C21  C22  1.3918  N21  110.7737   C23    0.0000
            C23  N21  1.3351  C22  105.1514   C21   -0.0000
            N22  C21  1.3815  C22  104.8582   N21    0.0000
            H21  C23  1.0872  N21  126.0493   C22 -180.0000
            H22  N22  1.0130  C21  126.1042   C22  180.0000
            H23  C22  1.0878  N21  121.5378   C23  180.0000
            H24  C21  1.0868  C22  132.6904   N21 -180.0000
            X21  N21  1.0000  C22   90.0000   C21  180.0000
            0 1
        '''

        atom_name = [
            'NE2', 'CD2', 'CG',  'CE1', 'ND1', 'HE1', 'HD1', 'HD2', 'HG'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def monomer_b_aromatic_zmatrix(self):

      zmatrix = '''\
            X21   :1  DISTANCE  :2   ANGLE  :3   90.0000
            C21  X21  1.1134    :1   90.0000   :2   0.0000
            C22  C21  1.3918  X21   60.0000    :1  DIHEDRAL
            N21  C22  1.3730  C21  110.7737   X21   0.0000
            C23  N21  1.3351  C22  105.1514   C21   -0.0000
            N22  C21  1.3815  C22  104.8582   N21    0.0000
            H21  C23  1.0872  N21  126.0493   C22 -180.0000
            H22  N22  1.0130  C21  126.1042   C22  180.0000
            H23  C22  1.0878  N21  121.5378   C23  180.0000
            H24  C21  1.0868  C22  132.6904   N21 -180.0000
            0 1
        '''

      atom_name = [
        'CG', 'CD2', 'NE2', 'CE1', 'ND1', 'HE1', 'HD1', 'HD2', 'HG'
      ]

      return textwrap.dedent(zmatrix), atom_name

