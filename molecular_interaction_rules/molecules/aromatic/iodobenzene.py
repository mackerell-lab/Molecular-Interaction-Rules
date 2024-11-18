#!/usr/bin/env python3
#
# Lennard-Jones-Drill-2: IodoBenzene
# ----------------------------------

# Imports
# -------
import textwrap


class IodoBenzene(object):

    def __init__(self):

        self.resi_name = 'IODB'

    def get_monomer_a_species(self):

      '''

      Get the Monomer A Species

      '''

      monomer_a_species = {
          'I1': self.monomer_a_halogen_zmatrix(),
          'RC1': self.monomer_a_centroid_zmatrix()
      }

      return monomer_a_species

    def get_monomer_b_species(self):

      monomer_b_species = {
        'I1': self.monomer_b_halogen_zmatrix(),
        'RC1': self.monomer_b_centroid_zmatrix()
      }

      return monomer_b_species

    def monomer_a_centroid_zmatrix(self):

      zmatrix = '''\
            X11
            C11  X11  1.3940
            C12  C11  1.3774 X11   60.0000
            C13  C12  1.3774 C11  120.0000 X11    0.0000
            C14  C13  1.3774 C12  120.0000 C11    0.0000
            C15  C14  1.3774 C13  120.0000 C12    0.0000
            C16  C15  1.3774 C14  120.0000 C13    0.0000
            I11  C11  2.2114 C12  120.0000 C13  180.0000
            H11  C12  1.0756 C11  120.0000 C13  180.0000
            H12  C13  1.0756 C12  120.0000 C11  180.0000
            H13  C14  1.0756 C13  120.0000 C12  180.0000
            H14  C15  1.0756 C14  120.0000 C13  180.0000
            H15  C16  1.0756 C15  120.0000 C11  180.0000
            0 1
        '''

      atom_name = [
        'C6', 'C5', 'C4', 'C3', 'C2', 'C1', 'I6', 'H5', 'H4', 'H3', 'H2', 'H1'
      ]

      return textwrap.dedent(zmatrix), atom_name


    def monomer_a_halogen_zmatrix(self):

      zmatrix = '''\
          I11
          C11  I11  2.2114
          C12  C11  1.3986   I11  120.0000
          C13  C12  1.3774   C11  120.0000   I11  -180.0000
          C14  C13  1.3774   C12  120.0000   C11     0.0000
          C15  C14  1.3774   C13  120.0000   C12     0.0000
          C16  C11  1.3774   C12  120.0000   C13    -0.0000
          H11  C13  1.0756   C12  120.0000   C11   180.0000
          H12  C14  1.0756   C13  120.0000   C12  -180.0000
          H13  C15  1.0756   C14  120.0000   C13   180.0000
          H14  C16  1.0756   C11  120.0000   I11    -0.0000
          H15  C12  1.0756   C11  120.0000   I11     0.0000
          X11  I11  1.0000   C11   90.0000   C12   180.0000
          0 1
       '''
      atom_name = [
        'I6', 'C6', 'C5', 'C4', 'C3', 'C2', 'C1', 'H5', 'H4', 'H3', 'H2', 'H1'
      ]

      return textwrap.dedent(zmatrix), atom_name

    def monomer_b_halogen_zmatrix(self):

      zmatrix = '''\
            I21  :1  DISTANCE   :2   ANGLE     :3   90.0000
            X21  I21  1.0000    :1  90.0000    :2    0.0000
            C21  I21  2.2114    X21  90.0000   :1   DIHEDRAL
            C22  C21  1.3986   I21  120.0000   :1     0.0000
            C23  C22  1.3774   C21  120.0000   I21  -180.0000
            C24  C23  1.3774   C22  120.0000   C21     0.0000
            C25  C24  1.3774   C23  120.0000   C22     0.0000
            C26  C21  1.3774   C22  120.0000   C23    -0.0000
            H21  C23  1.0756   C22  120.0000   C21   180.0000
            H22  C24  1.0756   C23  120.0000   C22  -180.0000
            H23  C25  1.0756   C24  120.0000   C23   180.0000
            H24  C26  1.0756   C21  120.0000   I21    -0.0000
            H25  C22  1.0756   C21  120.0000   I21     0.0000
            0 1
        '''

      atom_name = [
        'I6', 'C6', 'C5', 'C4', 'C3', 'C2', 'C1', 'H5', 'H4', 'H3', 'H2', 'H1'
      ],

      return textwrap.dedent(zmatrix), atom_name

    def monomer_b_centroid_zmatrix(self):

      zmatrix = '''\
            X21   :1  DISTANCE  :2   ANGLE      :3     90.0000
            C21  X21  1.3940    :1   90.0000    :2  180.0000
            C22  C21  1.3774 X21   60.0000      :1   90.0000
            C23  C22  1.3774 C21  120.0000 X21   DIHEDRAL
            C24  C23  1.3774 C22  120.0000 C21    0.0000
            C25  C24  1.3774 C23  120.0000 C22    0.0000
            C26  C25  1.3774 C24  120.0000 C23    0.0000
            I21  C21  2.2114 C22  120.0000 C23  180.0000
            H21  C22  1.0756 C21  120.0000 C23  180.0000
            H22  C23  1.0756 C22  120.0000 C21  180.0000
            H23  C24  1.0756 C23  120.0000 C22  180.0000
            H24  C25  1.0756 C24  120.0000 C23  180.0000
            H25  C26  1.0756 C25  120.0000 C21  180.0000
            0 1
        '''

      atom_name = [
        'C6', 'C5', 'C4', 'C3', 'C2', 'C1', 'I6', 'H5', 'H4', 'H3', 'H2', 'H1'
      ]

      return textwrap.dedent(zmatrix), atom_name

