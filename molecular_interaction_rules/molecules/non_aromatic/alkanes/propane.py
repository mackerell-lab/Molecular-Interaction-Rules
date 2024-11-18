# Imports
# -------

import textwrap

class Propane(object):

    def __init__(self):

        self.resi_name = 'PRPA'

    def get_monomer_a_species(self):

        '''

        Get the Monomer A Species

        '''

        monomer_a_species = {
            'H1': self.get_monomer_a_terminal_hydrogen_zmatrix(),
            'H2': self.get_monomer_a_center_hydrogen_zmatrix(),
        }

        return monomer_a_species

    def get_monomer_b_species(self):

        monomer_b_species = {
          'H1': self.get_monomer_b_terminal_hydrogen_zmatrix(),
          'H2': self.get_monomer_b_terminal_hydrogen_zmatrix()
        }

        return monomer_b_species

    def get_monomer_a_terminal_hydrogen_zmatrix(self):

        zmatrix = '''\
            H11
            C12 H11 1.1026
            C11 C12 1.5332 H11 110.7849
            H12 C12 1.1026 C11  110.7849   H11  119.4559
            H13 C12 1.1026 C11  110.7849   H11 -119.4559
            C13 C11 1.5332 C12  112.0949   H11   59.7279
            H14 C13 1.1026 C11  111.6593   C12   59.7279
            H15 C13 1.1026 C11  111.6593   C12  -59.7279
            H16 C13 1.1014 C11  111.6593   C12  178.5023
            H17 C11 1.1034 C12  109.5193   H11 -178.5023
            H18 C11 1.1034 C12  109.5193   H11  -62.0416
            X11 H11 1.0000 C12   90.0000   C11  180.0000
            0 1
        '''

        atom_name = [
          'H11', 'C1', 'C2', 'H12', 'H13', 'C3', 'H31', 'H32', 'H33', 'H21', 'H22'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_b_terminal_hydrogen_zmatrix(self):

        zmatrix = '''\
              H21   :1 DISTANCE  :2 ANGLE :3  DIHEDRAL
              X21 H21 1.0000  :1  90.0000   :2   0.0000
              C12 H11 1.1026 X21  90.0000   :1  180.0000
              C11 C12 1.5332 H11 110.7849   1  180.0000
              H12 C12 1.1026 C11  110.7849   H11  119.4559
              H13 C12 1.1026 C11  110.7849   H11 -119.4559
              C13 C11 1.5332 C12  112.0949   H11   59.7279
              H14 C13 1.1026 C11  111.6593   C12   59.7279
              H15 C13 1.1026 C11  111.6593   C12  -59.7279
              H16 C13 1.1014 C11  111.6593   C12  178.5023
              H17 C11 1.1034 C12  109.5193   H11 -178.5023
              H18 C11 1.1034 C12  109.5193   H11  -62.0416
              X11 H11 1.0000 C12   90.0000   C11  180.0000
              0 1
          '''

        atom_name = [
          'H11', 'C1', 'C2', 'H12', 'H13', 'C3', 'H31', 'H32', 'H33', 'H21', 'H22'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_a_center_hydrogen_zmatrix(self):

        zmatrix = '''\
            H11
            C11 H11 1.1034
            C12 C11 1.5332 H11 109.5193
            H12 C12 1.1026 C11 110.7849    H11 -178.5023
            H13 C12 1.1026 C11  110.7849   H12  119.4559
            H14 C12 1.1026 C11  110.7849   H12 -119.4559
            C13 C11 1.5332 C12  112.0949   H12   59.7279
            H15 C13 1.1026 C11  111.6593   C12   59.7279
            H16 C13 1.1026 C11  111.6593   C12  -59.7279
            H17 C13 1.1014 C11  111.6593   C12  178.5023
            H18 C11 1.1034 C12  109.5193   H12  -62.0416
            X11 H11 1.0000 C11   90.0000   C12  180.0000
            0 1
        '''

        atom_name = [
          'H21', 'C2', 'C1', 'H11', 'H12', 'H13', 'C3', 'H31', 'H32', 'H33', 'H22'
        ]

        return textwrap.dedent(zmatrix), atom_name


    def get_monomer_b_center_hydrogen_zmatrix(self):

      zmatrix = '''\
            H21  :1 DISTANCE  :2 ANGLE :3  DIHEDRAL
            X21 H21 1.0000  :1  90.0000   :2   0.0000
            C11 H11 1.1034  X21  90.0000   :1  180.0000
            C12 C11 1.5332 H11 109.5193    :1  180.0000
            H12 C12 1.1026 C11 110.7849    H11 -178.5023
            H13 C12 1.1026 C11  110.7849   H12  119.4559
            H14 C12 1.1026 C11  110.7849   H12 -119.4559
            C13 C11 1.5332 C12  112.0949   H12   59.7279
            H15 C13 1.1026 C11  111.6593   C12   59.7279
            H16 C13 1.1026 C11  111.6593   C12  -59.7279
            H17 C13 1.1014 C11  111.6593   C12  178.5023
            H18 C11 1.1034 C12  109.5193   H12  -62.0416
            X11 H11 1.0000 C11   90.0000   C12  180.0000
            0 1
        '''

      atom_name = [
        'H21', 'C2', 'C1', 'H11', 'H12', 'H13', 'C3', 'H31', 'H32', 'H33', 'H22'
      ]

      return textwrap.dedent(zmatrix), atom_name

