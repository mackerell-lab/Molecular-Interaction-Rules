#!/usr/bin/env python3
#
# Lennard-Jones-Drill-2: Isoxazole
# ---------------------------------

# Imports
# -------
import textwrap

class Isoxazole(object):

    def __init__(self):

        self.resi_name = 'isox'

    def get_monomer_a_species(self):

        monomer_a_species = {
            'RC1': self.monomer_a_aromatic_zmatrix(),
            'O1': self.monomer_a_oxygen_zmatrix(),
            'N1': self.get_monomer_a_nitrogen()
        }

        return monomer_a_species

    def get_monomer_b_species(self):

        monomer_b_species = {
            'RC1': self.monomer_b_pi_stack_zmatrix(),
        }

        return monomer_b_species

    def monomer_a_aromatic_zmatrix(self):

        zmatrix = '''\
          X11
          O11 X11 1.1000
          N11 O11 1.3936   X11   60.0000
          C11 N11 1.3368   O11  105.3771   X11  180.0000
          C12 C11 1.4240   N11  111.9773   O11   -0.0000
          C13 O11 1.3580   N11  109.3403   C11   -0.0000
          H11 C12 1.0866   C11  128.8980   N11 -180.0000
          H12 C13 1.0873   O11  115.7339   N11 -180.0000
          H13 C11 1.0888   N11  118.3782   O11 -180.0000
          0 1
        '''

        atom_name = [
            'O5', 'N4', 'C3', 'C2', 'C1', 'H2', 'H3', 'H1'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def monomer_a_oxygen_zmatrix(self):

        zmatrix = '''\
          O11
          N11 O11 1.3936
          C11 N11 1.3368   O11  105.3771
          C12 C11 1.4240   N11  111.9773   O11   -0.0000
          H11 C12 1.0866   C11  128.8980   N11 -180.0000
          C13 O11 1.3580   N11  109.3403   C11   -0.0000
          H12 C13 1.0873   O11  115.7339   N11 -180.0000
          H13 C11 1.0888   N11  118.3782   O11 -180.0000
          X11 O11  1.0000  N11   90.0000   C11  180.0000
          0 1
        '''

        atom_name = [
            'O5', 'N4', 'C3', 'C2', 'C1', 'H2', 'H3', 'H1'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_a_nitrogen(self):

        zmatrix = '''\
            O11
            N11 O11 1.3936
            C11 N11 1.3368   O11  105.3771
            C12 C11 1.4240   N11  111.9773   O11   -0.0000
            H11 C12 1.0866   C11  128.8980   N11 -180.0000
            C13 O11 1.3580   N11  109.3403   C11   -0.0000
            H12 C13 1.0873   O11  115.7339   N11 -180.0000
            H13 C11 1.0888   N11  118.3782   O11 -180.0000
            X11 O11  1.0000  N11   90.0000   C11  180.0000
            0 1
          '''

        atom_name = [
          'O5', 'N4', 'C3', 'C2', 'C1', 'H2', 'H3', 'H1'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def monomer_b_pi_stack_zmatrix(self):

        zmatrix = '''\
          X21  :1  DISTANCE  :2   ANGLE    :3   90.0000
          O21 X21 1.1000    :1   90.0000   :2  180.0000
          N21 O21 1.3936   X21   60.0000   :1   90.0000
          C21 N21 1.3368   O21  105.3771   X21  DIHEDRAL
          C22 C21 1.4240   N21  111.9773   O21   -0.0000
          C23 O21 1.3580   N21  109.3403   C21   -0.0000
          H21 C22 1.0866   C21  128.8980   N21 -180.0000
          H22 C23 1.0873   O21  115.7339   N21 -180.0000
          H23 C21 1.0888   N21  118.3782   O21 -180.0000
          0 1
        '''

        atom_name = [
          'O5', 'N4', 'C3', 'C2', 'C1', 'H2', 'H3', 'H1'
        ]

        return textwrap.dedent(zmatrix), atom_name
