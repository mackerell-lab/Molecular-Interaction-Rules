#!/usr/bin/env python3
#
# Lennard-Jones-Drill-2: Imidazolium
# ----------------------------------

# Imports
# -------
import textwrap

class Imidazolium(object):

    def __init__(self):

        self.resi_name = 'imim'

    def get_monomer_a_species(self):

        '''

        Get the Monomer A Species

        '''

        monomer_a_species = {
            'RC1': self.monomer_a_aromatic_zmatrix(),
            'H1': self.monomer_a_hydrogen_zmatrix()
        }

        return monomer_a_species

    def get_monomer_b_species(self):

        monomer_b_species = {
            'RC1': self.monomer_b_aromatic_zmatrix(),
            'H1': self.monomer_b_hydrogen_zmatrix()
        }

        return monomer_b_species

    def monomer_a_aromatic_zmatrix(self):

        zmatrix = '''\
            X11
            C11  X11  1.1134
            C12  C11  1.3918  X11   59.0000
            N11  C12  1.3834  C11  110.7737   X11    0.0000
            C13  N11  1.3351  C12  105.1514   C11   -0.0000
            N12  C11  1.3815  C12  104.8582   N11    0.0000
            H11  C13  1.0872  N11  126.0493   C12 -180.0000
            H12  N12  1.0130  C11  126.1042   C12  180.0000
            H13  C12  1.0878  N11  121.5378   C13  180.0000
            H14  C11  1.0868  C12  132.6904   N11 -180.0000
            H15  N11  1.0130  C12  126.1042   C11  180.0000
            1 1
        '''

        atom_name = [
            'CG', 'CD2', 'NE2', 'CE1', 'ND1', 'HE1', 'HD1', 'HD2', 'HG', 'HE2'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def monomer_b_aromatic_zmatrix(self):

      zmatrix = '''\
            X21    :1  DISTANCE  :2   ANGLE  :3   90.0000
            C21  X21  1.1134     :1   90.0000   :2  180.0000
            C22  C21  1.3918  X21   59.0000     :1   DIHEDRAL
            N21  C22  1.3834  C21  110.7737   X21    0.0000
            C23  N21  1.3351  C22  105.1514   C21   -0.0000
            N22  C21  1.3815  C22  104.8582   N21    0.0000
            H21  C23  1.0872  N21  126.0493   C22 -180.0000
            H22  N22  1.0130  C21  126.1042   C22  180.0000
            H23  C22  1.0878  N21  121.5378   C23  180.0000
            H24  C21  1.0868  C22  132.6904   N21 -180.0000
            H25  N21  1.0130  C22  126.1042   C21  180.0000
            1 1
        '''

      atom_name = [
        'CG', 'CD2', 'NE2', 'CE1', 'ND1', 'HE1', 'HD1', 'HD2', 'HG', 'HE2'
      ]

      return textwrap.dedent(zmatrix), atom_name

    def monomer_a_nitrogen_zmatrix(self):

        zmatrix = '''\
            N11
            C12  N11  1.3834
            C11  C12  1.3918  N11  110.7737
            C13  N11  1.3351  C12  105.1514   C11   -0.0000
            N12  C11  1.3815  C12  104.8582   N11    0.0000
            H11  C13  1.0872  N11  126.0493   C12 -180.0000
            H12  N12  1.0130  C11  126.1042   C12  180.0000
            H13  C12  1.0878  N11  121.5378   C13  180.0000
            H14  C11  1.0868  C12  132.6904   N11 -180.0000
            H15  N11  1.0130  C12  126.1042   C11  180.0000
            X11  N11  1.0000  C12   90.0000   C11  180.0000
            1 1
        '''

        atom_name = [
          'NE2', 'CD2', 'CG',  'CE1', 'ND1', 'HE1', 'HD1', 'HD2', 'HG', 'HE2'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def monomer_a_hydrogen_zmatrix(self):

        zmatrix = '''\
            H11
            C13  H11  1.0782
            N11  C13  1.3351  H11  126.1042
            C12  N11  1.3834  C13  105.1514   H11 -180.0000
            C11  C12  1.3918  N11  110.7737   C13    0.0000
            N12  C11  1.3815  C12  104.8582   N11    0.0000
            H12  N12  1.0130  C11  126.1042   C12  180.0000
            H13  C12  1.0872  N11  121.5378   C13  180.0000
            H14  C11  1.0872  C12  132.6904   N11 -180.0000
            H15  N11  1.0130  C12  126.1042   C11  180.0000
            X11  N11  1.0000  C12   90.0000   C11  180.0000
            1 1
        '''

        atom_name = [
          'HE1', 'CE1', 'NE2', 'CD2', 'CG', 'ND1','HD1', 'HD2', 'HG', 'HE2'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def monomer_b_hydrogen_zmatrix(self):

        zmatrix = '''\
              H21  :1  DISTANCE  :2 ANGLE      :3   90.0000
              X21  H21 1.0000    :1 90.0000      :2    0.0000
              C23  H21  1.0782  X21  90.0000     :1  DIHEDRAL
              N21  C23  1.3351  H21  126.1042    :1    0.0000
              C22  N21  1.3834  C23  105.1514   H21 -180.0000
              C21  C22  1.3918  N21  110.7737   C23    0.0000
              N22  C21  1.3815  C22  104.8582   N21    0.0000
              H22  N22  1.0130  C21  126.1042   C22  180.0000
              H23  C22  1.0872  N21  121.5378   C23  180.0000
              H24  C21  1.0872  C22  132.6904   N21 -180.0000
              H25  N21  1.0130  C22  126.1042   C21  180.0000
              1 1
          '''

        atom_name = [
          'HE1', 'CE1', 'NE2', 'CD2', 'CG', 'ND1','HD1', 'HD2', 'HG', 'HE2'
        ]

        return textwrap.dedent(zmatrix), atom_name


