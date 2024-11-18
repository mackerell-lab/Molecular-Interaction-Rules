# Imports
# -------
import textwrap

class Methanol(object):

    def __init__(self):

        self.resi_name = 'MEOH'

    def get_monomer_a_species(self):

        '''

        Get the Monomer A Z-Matrices

        '''

        monomer_a_species = {
            'O1': self.get_oxygen_zmatrix(),
            'H1': self.hydrogen_acceptor_zmatrix(),
        }

        return monomer_a_species

    def get_monomer_b_species(self):

        '''

        Get the Monomer A Z-Matrices

        '''

        monomer_b_species = {
            'O1': self.oxygen_acceptor_b_zmatrix(),
            'H1': self.monomer_b_hydrogen_acceptor_zmatrix(),
        }

        return monomer_b_species

    def get_oxygen_zmatrix(self):

        zmatrix = '''\
            O11
            H11  O11  0.9657
            C11  O11  1.4349  H11  126.2746
            H12  C11  1.1029  O11  111.8699  H11    0.0000
            H13  C11  1.1029  O11  111.8699  H11  122.9683
            H14  C11  1.1029  O11  111.8699  H11 -118.5158
            X11  O11  1.0000  H11   90.0000  C11  180.0000
            0 1
        '''

        atom_name = [
          'OG', 'HG1', 'CB', 'HB1', 'HB2', 'HB3'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def hydrogen_acceptor_zmatrix(self):

        zmatrix = '''\
            H11
            O11  H11  0.9657
            C11  O11  1.4349  H11  107.5890
            H12  C11  1.1029  O11  111.8699  H11    0.0000
            H13  C11  1.1029  O11  111.8699  H11  122.9683
            H14  C11  1.1029  O11  111.8699  H11 -118.5158
            X11  H11  1.0000  O11   90.0000  C11  180.0000
            0 1
        '''

        atom_name = [
          'HG1', 'OG', 'CB', 'HB1', 'HB2', 'HB3'
        ]

        return textwrap.dedent(zmatrix), atom_name


    def oxygen_acceptor_b_zmatrix(self):

        zmatrix = '''\
            O21   :1  DISTANCE  :2 ANGLE      :3    DIHEDRAL
            H21  O21  0.9657  :1  236.0000    :2    0.0000
            X21  H21  1.0000   O21   90.0000  :1    90.0000
            C21  O21  1.4349   H21 111.8699    X21   0.0000
            H22  C21  1.1029  O21  111.8699  H21    0.0000
            H23  C21  1.1029  O21  111.8699  H21  122.9683
            H24  C21  1.1029  O21  111.8699  H21 -118.5158
            0 1
        '''

        atom_name = [
          'OG', 'HG1', 'CB', 'HB2', 'HB2', 'HB3'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def monomer_b_hydrogen_acceptor_zmatrix(self):

      zmatrix = '''\
            H21   :1  DISTANCE  :2   ANGLE     :3   90.0000
            X21  H21    1.0000  :1   90.0000   :2   180.0000
            O21  H21    0.9657  X21  90.0000   :1   180.0000
            C21  O21    1.4349  H21 120.0000  X21  DIHEDRAL
            H22  C21    1.1029  O21   111.8699  H21    0.0000
            H23  C21    1.1029  O21   111.8699  H21  122.9683
            H24  C21    1.1029  O21   111.8699  H21 -118.5158
            0 1
        '''

      atom_name = [
        'HG1', 'OG', 'CB', 'HB1', 'HB2', 'HB3'
      ]

      return textwrap.dedent(zmatrix), atom_name
