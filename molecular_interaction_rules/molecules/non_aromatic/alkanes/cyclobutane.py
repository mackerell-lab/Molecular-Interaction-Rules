# Imports
# -------

import textwrap

class Cyclobutane(object):

    def __init__(self):

        self.resi_name = 'CBU'

    def get_monomer_a_species(self):

        '''

        Get the Monomer A Species

        '''

        monomer_a_species = {
          'H1': self.monomer_a_hydrogen_zmatrix()
        }

        return monomer_a_species

    def get_monomer_b_species(self):

        monomer_b_species = {
            'H1': self.monomer_b_hydrogen_zmatrix()
        }

        return monomer_b_species

    def monomer_a_hydrogen_zmatrix(self):

        zmatrix = '''\
            H11
            C11  H11  1.1001
            C12  C11  1.5573  H11  118.5122
            C13  C12  1.5573  C11   87.8052  H11 -143.6662
            C14  C11  1.5573  C12   87.8052  C13  -22.1452
            H12  C14  1.1015  C11  110.6146  C12  -89.0794
            H13  C14  1.1001  C11  118.5122  C12  143.6662
            H14  C13  1.1015  C12  110.6146  C11  -89.0794
            H15  C13  1.1001  C12  118.5122  C11  143.6662
            H16  C12  1.1001  C11  118.5122  C14 -143.6662
            H17  C12  1.1015  C11  110.6145  C14   89.0793
            H18  C11  1.1015  C12  110.6146  C13   89.0794
            X11  H11  1.0000  C11   90.0000  C12  180.0000
            0 1
        '''

        atom_name = [
          'H11', 'C1', 'C2', 'C3', 'C4', 'H41', 'H42', 'H31', 'H32', 'H21', 'H22', 'H12'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def monomer_b_hydrogen_zmatrix(self):

      zmatrix = '''\
            H21  :1  DISTANCE  :2 ANGLE :3  DIHEDRAL
            X21 H21 1.0000  :1  90.0000  :2  0.0000
            C21 H21 1.1013 X21  90.0000  :2  180.0000
            C22 C21 1.5345 H21 111.6112  :1  180.0000
            C23  C22  1.5573  C21   87.8052  H21 -143.6662
            C24  C21  1.5573  C22   87.8052  C23  -22.1452
            H22  C24  1.1015  C21  110.6146  C22  -89.0794
            H23  C24  1.1001  C21  118.5122  C22  143.6662
            H24  C23  1.1015  C22  110.6146  C21  -89.0794
            H25  C23  1.1001  C22  118.5122  C21  143.6662
            H26  C22  1.1001  C21  118.5122  C24 -143.6662
            H27  C22  1.1015  C21  110.6145  C24   89.0793
            H28  C21  1.1015  C22  110.6146  C23   89.0794
            0 1
        '''

      atom_name = [
        'H11', 'C1', 'C2', 'C3', 'C4', 'H41', 'H42', 'H31', 'H32', 'H21', 'H22', 'H12'
      ]

      return textwrap.dedent(zmatrix), atom_name
