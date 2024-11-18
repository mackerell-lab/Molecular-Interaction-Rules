# Imports
# -------

import textwrap

class Methoxyethene(object):

    def __init__(self):

        self.resi_name = 'MOET'

    def get_monomer_a_species(self):

        monomer_a_species = {
            'O1': self.get_monomer_a_oxygen_zmatrix(),
            'H1': self.get_monomer_a_hydrogen_zmatrix(),
            'H2': self.get_monomer_a_hydrogen_inner_zmatrix()
        }

        return monomer_a_species

    def get_monomer_b_species(self):

      '''

      Get the Monomer A Z-Matrices

      '''

      monomer_b_species = {
          'O1': self.get_monomer_b_oxygen_zmatrix(),
          'H1': self.get_monomer_b_hydrogen_zmatrix(),
          'H2': self.get_monomer_b_hydrogen_inner_zmatrix()
      }

      return monomer_b_species

    def get_monomer_a_hydrogen_inner_zmatrix(self):

        zmatrix = '''\
              H12
              C11 H12 1.0912
              C12 C11 1.3486 H12 121.3771
              H11 C12 1.0912 C11 121.3771 H12    0.0000
              O11 C11 1.3716 C12 122.1594 H11 -180.0000
              C13 O11 1.3716 C11 113.4660 C12 -169.9762
              H14 C13 1.0966 O11 110.7661 C11   63.2058
              H15 C13 1.0966 O11 110.7661 C11  -58.7758
              H16 C13 1.0966 O11 110.7661 C11 -177.7054
              X11 H11 1.0000 C12  90.0000 C11  180.0000
              0 1
        '''

        atom_name = [
          'H21', 'C2', 'C1', 'H11', 'O3', 'C4','H12','H41', 'H42', 'H43'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_a_oxygen_zmatrix(self):

        zmatrix = '''\
            O11
            C11 O11 1.3716
            C12 C11 1.3486 O11 122.1594
            C13 O11 1.3716 C11 113.4660 C12 -169.9762
            H11 C12 1.0912 C11 121.3771 O11 -180.0000
            H12 C12 1.0912 C11 121.3771 O11    0.0000
            H14 C13 1.0966 O11 110.7661 C11   63.2058
            H15 C13 1.0966 O11 110.7661 C11  -58.7758
            H16 C13 1.0966 O11 110.7661 C11 -177.7054
            X11 O11 1.0000 C11  90.0000 C12  180.0000
            0 1
        '''

        atom_name = [
          'O3', 'C2', 'C1', 'C4', 'H11', 'H12', 'H21', 'H41', 'H42', 'H43'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_b_oxygen_zmatrix(self):

      zmatrix = '''\
            O21  :1 DISTANCE :2 ANGLE    :3 90.0000
            C21 O21 1.4635  :1 126.2746  :2  DIHEDRAL
            C22 C21 1.5656 O21 107.5604  :1  0.0000
            C23 O21 1.3716 C21 113.4660 C22 -169.9762
            H21 C22 1.0912 C21 121.3771 O21 -180.0000
            H22 C22 1.0912 C21 121.3771 O21    0.0000
            H23 C21 1.0912 C22 121.3771 H21    0.0000
            H24 C23 1.0966 O21 110.7661 C21   63.2058
            H25 C23 1.0966 O21 110.7661 C21  -58.7758
            H26 C23 1.0966 O21 110.7661 C21 -177.7054
            X21 O21 1.0000 C21  90.0000 C22  180.0000
            0 1
        '''

      atom_name = [
        'O3', 'C2', 'C1', 'C4', 'H11', 'H12', 'H21', 'H41', 'H42', 'H43'
      ]

      return textwrap.dedent(zmatrix), atom_name

    def get_monomer_a_hydrogen_zmatrix(self):

        zmatrix = '''\
            H11
            C12 H11 1.0912
            C11 C12 1.3486 H11 121.3771
            O11 C11 1.3716 C12 122.1594 H11 -180.0000
            C13 O11 1.3716 C11 113.4660 C12 -169.9762
            H12 C12 1.0912 C11 121.3771 O11    0.0000
            H14 C13 1.0966 O11 110.7661 C11   63.2058
            H15 C13 1.0966 O11 110.7661 C11  -58.7758
            H16 C13 1.0966 O11 110.7661 C11 -177.7054
            X11 H11 1.0000 C12  90.0000 C11  180.0000
            0 1
        '''

        atom_name = [
          'H11', 'C1', 'C2', 'O3', 'C4','H12', 'H21', 'H41', 'H42', 'H43'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_b_hydrogen_zmatrix(self):

        zmatrix = '''\
              H21  :1  DISTANCE :2   ANGLE   :3 DIHEDRAL
              X21 H21  1.0000   :1  90.0000  :2    0.0000
              C22 H21 1.0912  X21  90.0000     :2  180.0000
              C21 C22 1.3486 H21 121.3771      :1  180.0000
              O21 C21 1.3716 C22 122.1594 H21 -180.0000
              C23 O21 1.3716 C21 113.4660 C22 -169.9762
              H22 C22 1.0912 C21 121.3771 O21    0.0000
              H24 C23 1.0966 O21 110.7661 C21   63.2058
              H25 C23 1.0966 O21 110.7661 C21  -58.7758
              H26 C23 1.0966 O21 110.7661 C21 -177.7054
              X21 H21 1.0000 C22  90.0000 C21  180.0000
              0 1
          '''

        atom_name = [
          'H11', 'C1', 'C2', 'O3', 'C4','H12', 'H21', 'H41', 'H42', 'H43'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_b_hydrogen_inner_zmatrix(self):

        zmatrix = '''\
            H22   :1  DISTANCE :2 ANGLE  :3 DIHEDRAL
            X21 H21   1.0000   :1  90.0000  :2    0.0000
            C21 H22 1.0912  X21  90.0000    :2  180.0000
            C22 C21 1.3486 H22 121.3771     :1  180.0000
            H21 C22 1.0912 C21 121.3771 H22    0.0000
            O21 C21 1.3716 C22 122.1594 H21 -180.0000
            C23 O21 1.3716 C21 113.4660 C22 -169.9762
            H24 C23 1.0966 O21 110.7661 C21   63.2058
            H25 C23 1.0966 O21 110.7661 C21  -58.7758
            H26 C23 1.0966 O21 110.7661 C21 -177.7054
            0 1
        '''

        atom_name = [
          'H21', 'C2', 'C1', 'H11', 'O3', 'C4','H12','H41', 'H42', 'H43'
        ]

        return textwrap.dedent(zmatrix), atom_name
