#!/usr/bin/env python3
#
# Lennard-Jones-Drill-2: Pyridine
# -------------------------------

# Imports
# -------
import textwrap

class Pyridine(object):

    def __init__(self):

        self.resi_name = 'PYR1'

    def get_monomer_a_species(self):

        '''

        Get the Monomer A Species

        '''

        monomer_a_species = {
            'RC1': self.monomer_a_aromatic_zmatrix()
        }

        return monomer_a_species

    def get_monomer_b_species(self):

        monomer_b_species = {
          'RC1': self.get_monomer_b_aromatic_zmatrix(),
        }

        return monomer_b_species

    def monomer_a_aromatic_zmatrix(self):

        zmatrix = '''\
            X11
            C11  X11  1.3940
            C12  C11  1.4045 X11   60.0000
            C13  C12  1.4045 C11  120.0000 X11    0.0000
            N11  C13  1.3208 C12  124.0000 C11    0.0000
            C14  N11  1.3208 C13  116.0000 C12    0.0000
            C15  C14  1.4045 N11  120.0000 C13    0.0000
            H11  C11  1.0756 C12  120.0000 C13  180.0000
            H12  C12  1.0756 C11  120.0000 C13  180.0000
            H13  C13  1.0756 C12  120.0000 C11  180.0000
            H14  C14  1.0756 N11  120.0000 C13  180.0000
            H15  C15  1.0756 C14  120.0000 C11  180.0000
            0 1
        '''

        atom_name = [
            'CG', 'CE1', 'CD1', 'NZ', 'CE2', 'CD2',
            'HG', 'HE1', 'HD1', 'HD2', 'HE2',
        ]

        return textwrap.dedent(zmatrix), atom_name

    def monomer_a_nitrogen_zmatrix(self):

        zmatrix = '''\
           N11
           C11  N11  1.3208
           C12  C11  1.4045  N11  124.0000
           C13  C12  1.4045  C11  120.0000  N11   0.0000
           C14  C13  1.4045  C12  120.0000  C11   0.0000
           C15  N11  1.3208  C11  116.0000  C12   0.0000
           H11  C11  1.0756  C12  120.0000  C13 180.0000
           H12  C12  1.0756  C13  120.0000  C14 180.0000
           H13  C13  1.0756  C14  120.0000  C15 180.0000
           H14  C14  1.0756  C15  120.0000  N11 180.0000
           H15  C15  1.0756  N11  120.0000  C11 180.0000
           X11  N11  1.0000  C11   90.0000  C12 180.0000
           0 1
        '''

        atom_name = [
            'NZ', 'CG', 'CE1', 'CD1', 'CE2', 'CD2',
            'HG', 'HE1', 'HD1', 'HD2', 'HE2',
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_b_aromatic_zmatrix(self):

        zmatrix = '''\
            X21   :1  DISTANCE  :2  ANGLE  :3   90.00000
            C21  X21  1.3940  :1  90.0000   :2  180.0000
            C22  C21  1.4045 X21   60.0000   :1  90.0000
            C23  C22  1.4045 C21  120.0000 X21   DIHEDRAL
            N21  C23  1.3208 C22  124.0000 C21    0.0000
            C24  N21  1.3208 C23  116.0000 C22    0.0000
            C25  C24  1.4045 N21  120.0000 C23    0.0000
            H21  C21  1.0756 C22  120.0000 C23  180.0000
            H22  C22  1.0756 C21  120.0000 C23  180.0000
            H23  C23  1.0756 C22  120.0000 C21  180.0000
            H24  C24  1.0756 N21  120.0000 C23  180.0000
            H25  C25  1.0756 C24  120.0000 C21  180.0000
            0 1
          '''

        atom_name = [
          'CG', 'CE1', 'CD1', 'NZ', 'CE2', 'CD2',
          'HG', 'HE1', 'HD1', 'HD2', 'HE2',
        ]

        return textwrap.dedent(zmatrix), atom_name



