#!/usr/bin/env python3
#
# Lennard-Jones-Drill-2: Tetrahydrofuran
# --------------------------------------

# Imports
# -------
import textwrap

class Tetrahydrofuran(object):

    def __init__(self):

        self.resi_name = 'THF'

    def get_monomer_a_species(self):

        '''

        Get the Monomer A Species

        '''

        monomer_a_species = {
            'O1': self.get_monomer_a_oxygen_zmatrix(),
            'H1': self.get_monomer_a_hydrogen_zmatrix()
        }

        return monomer_a_species

    def get_monomer_a_oxygen_zmatrix(self):

        zmatrix = '''\
          O11
          C11  O11  1.4476
          C12  C11  1.5321  O11  106.0969
          C13  C12  1.5376  C11  101.1722  O11   33.0711
          H11  C13  1.1020  C12  110.1183  C11   77.3973
          H12  C13  1.1001  C12  113.6317  C11 -160.2817
          C14  O11  1.4476  C11  109.1641  C12  -12.9455
          H13  C14  1.1004  O11  108.0713  C11 -135.1527
          H14  C14  1.1052  O11  108.9385  C11  106.6650
          H15  C12  1.1001  C11  112.6033  O11  154.7274
          H16  C12  1.1020  C11  110.3595  O11  -83.4967
          H17  C11  1.1052  C12  111.0412  C13  -85.1543
          H18  C11  1.1004  C12  113.6297  C13  151.6691
          X11  O11  1.0000  C11   90.0000  C12  180.0000
          0 1
        '''

        atom_name = [
          "O1'", "C1'", "C2'", "C3'", "H31'", "H32'", "C4'", "H41'", "H42'", "H21'", "H22'", "H11'", "H12'"
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_a_hydrogen_zmatrix(self):

        zmatrix = '''\
          H11
          C11  H11  1.1052
          O11  C11  1.4476  H11  112.6033
          C12  C11  1.5321  O11  106.0969  H11  106.6650
          C13  C12  1.5376  C11  101.1722  O11   33.0711
          H12  C13  1.1052  C12  110.1183  C11   77.3973
          H13  C13  1.1052  C12  113.6317  C11 -160.2817
          C14  O11  1.4476  C11  109.1641  C12  -12.9455
          H14  C14  1.1004  O11  108.0713  C11 -135.1527
          H15  C14  1.1052  O11  108.9385  C11  106.6650
          H16  C12  1.1052  C11  112.6033  O11  154.7274
          H17  C12  1.1052  C11  110.3595  O11  -83.4967
          H18  C11  1.1052  C12  111.0412  C13  -85.1543
          X11  H11  1.0000  C11   90.0000  O11  180.0000
          0 1
        '''

        atom_name = [
          "H12'", "C1'", "O1'", "C2'", "C3'", "H31'", "H32'", "C4'", "H41'", "H42'", "H21'", "H22'", "H11'",
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_b_species(self):

        monomer_b_species = {
          'O1': self.get_monomer_b_oxygen_zmatrix(),
          'H1': self.get_monomer_b_hydrogen_zmatrix()
        }

        return monomer_b_species

    def get_monomer_b_hydrogen_zmatrix(self):

        zmatrix = '''\
            H21   :1    DISTANCE    :2 ANGLE       :3 DIHEDRAL
            X21   H21   1.0000      :1  90.0000       :2    0.0000
            C21  H21  1.1052       X21  90.0000       :2  180.0000
            O21  C21  1.4476  H21  112.6033           :1  180.0000
            C22  C21  1.5321  O21  106.0969  H21  106.6650
            C23  C22  1.5376  C21  101.1722  O21   33.0711
            H22  C23  1.1052  C22  110.1183  C21   77.3973
            H23  C23  1.1052  C22  113.6317  C21 -160.2817
            C24  O21  1.4476  C21  109.1641  C22  -12.9455
            H24  C24  1.1004  O21  108.0713  C21 -135.1527
            H25  C24  1.1052  O21  108.9385  C21  106.6650
            H26  C22  1.1052  C21  112.6033  O21  154.7274
            H27  C22  1.1052  C21  110.3595  O21  -83.4967
            H28  C21  1.1052  C22  111.0412  C23  -85.1543
            0 1
          '''

        atom_name = [
          "H12'", "C1'", "O1'", "C2'", "C3'", "H31'", "H32'", "C4'", "H41'", "H42'", "H21'", "H22'", "H11'",
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_b_oxygen_zmatrix(self):

        zmatrix = '''\
            O21  :1  DISTANCE  :2   ANGLE    :3  DIHEDRAL
            C21  O21  1.4476   :1  126.9550  :2   0.00000
            C22  C21  1.5321   :1  126.9550  :2  180.0000
            C23  C22  1.5376  C21  101.1722  O21   33.0711
            H21  C23  1.1020  C22  110.1183  C21   77.3973
            H22  C23  1.1001  C22  113.6317  C21 -160.2817
            C24  O21  1.4476  C21  109.1641  C22  -12.9455
            H23  C24  1.1004  O21  108.0713  C21 -135.1527
            H24  C24  1.1052  O21  108.9385  C21  106.6650
            H25  C22  1.1001  C21  112.6033  O21  154.7274
            H26  C22  1.1020  C21  110.3595  O21  -83.4967
            H27  C21  1.1052  C22  111.0412  C23  -85.1543
            H28  C21  1.1004  C22  113.6297  C23  151.6691
            0 1
          '''

        atom_name = [
          "O1'", "C1'", "C2'", "C3'", "H31'", "H32'", "C4'", "H41'", "H42'", "H21'", "H22'", "H11'", "H12'"
        ]

        return textwrap.dedent(zmatrix), atom_name

