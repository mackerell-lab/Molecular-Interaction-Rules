#!/usr/bin/env python3
#
# Lennard-Jones-Drill-2: Acetone
# ------------------------------

# Imports
# -------

import textwrap

class Acetone(object):

    def __init__(self):

        self.resi_name = 'ACO '

    def get_monomer_a_species(self):

        '''

        Get the Monomer A Species

        '''

        monomer_a_species = {
            'O1': self.get_monomer_a_oxygen_donor_zmatrix(),
            'H1': self.get_monomer_a_hydrogen_carbon_zmatrix(),
            'C1': self.get_monomer_a_carbon_zmatrix()
        }

        return monomer_a_species

    def get_monomer_b_species(self):

        monomer_b_species = {
          'O1': self.get_monomer_b_carbonyl_oxygen()
        }

        return monomer_b_species

    def get_monomer_a_hydrogen_carbon_zmatrix(self):

        zmatrix = '''\
          H11
          C12  H11  1.1025
          C11  C12  1.5180  H11  109.8660
          C13  C11  1.5180  C12  116.5720  H11  180.0000
          O11  C11  1.2307  C12  121.7139  C13  180.0000
          H12  C12  1.1025  C11  109.8860  C13   58.8163
          H13  C12  1.1025  C11  109.8860  C13  -58.8163
          H14  C13  1.1025  C11  109.8860  C12   58.8163
          H15  C13  1.1025  C11  109.8860  C12  -58.8163
          H16  C13  1.1025  C11  109.8860  C12  180.0000
          X11  H11  1.0000  C12   90.0000  C11  180.0000
          0 1
        '''

        atom_name = [
          'H21', 'C2', 'C1', 'C3', 'O1', 'H21', 'H22', 'H23', 'H31', 'H32', 'H33',
        ]

        return textwrap.dedent(zmatrix), atom_name
    def get_monomer_a_oxygen_donor_zmatrix(self):

        zmatrix = '''\
            O11
            C11  O11  1.2307
            C12  C11  1.5180  O11  121.7139
            C13  C11  1.5180  C12  116.5719  O11  180.0000
            H11  C12  1.1025  C11  109.8860  O11  121.1836
            H12  C12  1.1025  C11  109.8860  O11 -121.1836
            H13  C12  1.1025  C11  109.8860  O11    0.0000
            H14  C13  1.1025  C11  109.8860  O11 -121.1836
            H15  C13  1.1025  C11  109.8860  O11  121.1836
            H16  C13  1.1025  C11  109.8860  O11    0.0000
            X11  O11  1.0000  C11   90.0000  C12  180.0000
            0 1
        '''

        atom_name = [
          'O1', 'C1', 'C2', 'C3', 'H21', 'H22', 'H23', 'H31', 'H32', 'H33',
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_a_carbon_zmatrix(self):

        zmatrix = '''\
            C11
            C12  C11  1.5180
            C13  C11  1.5180  C12  116.5720
            O11  C11  1.2307  C12  121.7139  C13  180.0000
            H11  C12  1.1025  C11  109.8860  C13  180.0000
            H12  C12  1.1025  C11  109.8860  C13   58.8163
            H13  C12  1.1025  C11  109.8860  C13  -58.8163
            H14  C13  1.1025  C11  109.8860  C12   58.8163
            H15  C13  1.1025  C11  109.8860  C12  -58.8163
            H16  C13  1.1025  C11  109.8860  C12  180.0000
            X11  C11  1.0000  C12   90.0000  C13  180.0000
            0 1
        '''

        atom_name = [
          'C1', 'O1', 'C2', 'C3', 'H21', 'H22', 'H23', 'H31', 'H32', 'H33',
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_b_carbonyl_oxygen(self):

        zmatrix = '''\
            O21  :1  DISTANCE   :2  ANGLE    :3    DIHEDRAL
            C21  O21  1.2307    :1  180.0000    :2   180.0000
            C22  C21  1.3986   O21  120.0000    :1     0.0000
            C23  C21  1.5180  C22  116.5719  O21  180.0000
            H21  C22  1.1025  C21  109.8860  O21  121.1836
            H22  C22  1.1025  C21  109.8860  O21 -121.1836
            H23  C22  1.1025  C21  109.8860  O21    0.0000
            H24  C23  1.1025  C21  109.8860  O21 -121.1836
            H25  C23  1.1025  C21  109.8860  O21  121.1836
            H26  C23  1.1025  C21  109.8860  O21    0.0000
            0 1
        '''

        atom_name = [
          'O1', 'C1', 'C2', 'C3', 'H21', 'H22', 'H23', 'H31', 'H32', 'H33',
        ]

        return textwrap.dedent(zmatrix), atom_name


