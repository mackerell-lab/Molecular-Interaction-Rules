#!/usr/bin/env python3

# Imports
# -------
import textwrap

class Benzene(object):

    def __init__(self):

        self.resi_name = 'BENZ'

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
          'H1': self.monomer_b_t_shaped_zmatrix(),
        }

        return monomer_b_species

    def monomer_a_aromatic_zmatrix(self):

        '''

        if X11 is first monomer a
        then :1 C11 :2 C12 :C13

        '''

        zmatrix = '''\
            X11
            C11  X11  1.3940
            C12  C11  1.3774 X11   60.0000
            C13  C12  1.3774 C11  120.0000 X11    0.0000
            C14  C13  1.3774 C12  120.0000 C11    0.0000
            C15  C14  1.3774 C13  120.0000 C12    0.0000
            C16  C15  1.3774 C14  120.0000 C13    0.0000
            H11  C11  1.0756 C12  120.0000 C13  180.0000
            H12  C12  1.0756 C11  120.0000 C13  180.0000
            H13  C13  1.0756 C12  120.0000 C11  180.0000
            H14  C14  1.0756 C13  120.0000 C12  180.0000
            H15  C15  1.0756 C14  120.0000 C13  180.0000
            H16  C16  1.0756 C15  120.0000 C11  180.0000
            0 1
        '''

        atom_name = [
            'CZ', 'CE1', 'CD1', 'CG', 'CD2', 'CE2',
            'HZ', 'HE1', 'HD1', 'HG', 'HD2', 'HE2',
        ]

        return textwrap.dedent(zmatrix), atom_name

    def monomer_a_hydrogen_zmatrix(self):

        zmatrix = '''\
              H11
              C11  H11  1.0756
              C12  C11  1.3774 H11  120.0000
              C13  C12  1.3774 C11  120.0000 H11  180.0000
              C14  C13  1.3774 C12  120.0000 C11    0.0000
              C15  C14  1.3774 C13  120.0000 C12    0.0000
              C16  C15  1.3774 C14  120.0000 C13    0.0000
              H12  C12  1.0756 C11  120.0000 C13  180.0000
              H13  C13  1.0756 C12  120.0000 C11  180.0000
              H14  C14  1.0756 C13  120.0000 C12  180.0000
              H15  C15  1.0756 C14  120.0000 C13  180.0000
              H16  C16  1.0756 C15  120.0000 C11  180.0000
              X11  H11  1.0000 C11   90.0000 C12    0.0000
              0 1
          '''

        atom_name = [
          'HZ', 'CZ', 'CE1', 'CD1', 'CG', 'CD2', 'CE2',
          'HE1', 'HD1', 'HG', 'HD2', 'HE2',
        ]

        return textwrap.dedent(zmatrix), atom_name

    def monomer_b_t_shaped_zmatrix(self):

        zmatrix = '''\
            H21   :1  DISTANCE  :2   ANGLE  :3   90.0000
            X21  H21  1.0000  :1   90.0000  :2    0.0000
            C21  H21  1.0756 X21   90.0000  :1  DIHEDRAL
            C22  C21  1.3774 H21  120.0000  :1    0.0000
            C23  C22  1.3774 C21  120.0000 H21  180.0000
            C24  C23  1.3774 C22  120.0000 C21    0.0000
            C25  C24  1.3774 C23  120.0000 C22    0.0000
            C26  C25  1.3774 C24  120.0000 C23    0.0000
            H22  C22  1.0756 C21  120.0000 C23  180.0000
            H23  C23  1.0756 C22  120.0000 C21  180.0000
            H24  C24  1.0756 C23  120.0000 C22  180.0000
            H25  C25  1.0756 C24  120.0000 C23  180.0000
            H26  C26  1.0756 C25  120.0000 C21  180.0000
            0 1
        '''

        atom_name = [
           'HZ', 'CZ', 'CE1', 'CD1', 'CG', 'CD2', 'CE2',
           'HE1', 'HD1', 'HG', 'HD2', 'HE2',
        ],

        return textwrap.dedent(zmatrix), atom_name

    def monomer_b_aromatic_zmatrix(self):

        zmatrix = '''\
            X21  H11  DISTANCE  :2   ANGLE  :3   90.0000
            C21  X21  1.3940  :1   90.0000  :2  180.0000
            C22  C21  1.3774 X21   60.0000  :1   90.0000
            C23  C22  1.3774 C21  120.0000 X21  DIHEDRAL
            C24  C23  1.3774 C22  120.0000 C21    0.0000
            C25  C24  1.3774 C23  120.0000 C22    0.0000
            C26  C25  1.3774 C24  120.0000 C23    0.0000
            H21  C21  1.0756 C22  120.0000 C23  180.0000
            H22  C22  1.0756 C21  120.0000 C23  180.0000
            H23  C23  1.0756 C22  120.0000 C21  180.0000
            H24  C24  1.0756 C23  120.0000 C22  180.0000
            H25  C25  1.0756 C24  120.0000 C23  180.0000
            H26  C26  1.0756 C25  120.0000 C21  180.0000
            0 1
        '''

        atom_name = [
          'CZ', 'CE1', 'CD1', 'CG', 'CD2', 'CE2',
          'HZ', 'HE1', 'HD1', 'HG', 'HD2', 'HE2',
        ],

        return textwrap.dedent(zmatrix), atom_name
