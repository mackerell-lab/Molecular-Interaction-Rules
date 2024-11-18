#!/usr/bin/env python3
#
# Lennard-Jones-Drill-2: Acetaldehyde
# -----------------------------------

# Imports
# -------

import textwrap

class Acetaldehyde(object):

    def __init__(self):

        self.resi_name = ''

    def get_monomer_a_species(self):

        '''

        Get the Monomer A Species

        '''

        monomer_a_species = {
            'O1': self.get_monomer_a_oxygen_donor_zmatrix(),
            'C1': self.get_carbonyl_electron_acceptor()
        }

        return monomer_a_species

    def get_monomer_b_species(self):

        monomer_b_species = {
        }

        return monomer_b_species

    def get_monomer_a_oxygen_donor_zmatrix(self):

        zmatrix = '''\
            O11
            C11 O11 1.2261
            C12 C11 1.5077 O11 124.4022
            H11 C12 1.1028 C11 109.2276 O11    0.0000
            H12 C12 1.1028 C11 109.2276 O11  121.6416
            H13 C12 1.1028 C11 109.2276 O11 -121.6416
            H14 C11 1.1025 C12 115.6994 H11  180.0000
            X11 O11 1.0000 C11  90.0000 C12  180.0000
            0 1
        '''

        atom_name = [
            'O', 'C', 'CB', 'HB1', 'HB2', 'HB3', 'HA'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_carbonyl_electron_acceptor(self):

        zmatrix = '''\
            C11
            O11 C11 1.2261
            C12 C11 1.5077 O11 124.4022
            H11 C12 1.1028 C11 109.2276 O11    0.0000
            H12 C12 1.1028 C11 109.2276 O11  121.6416
            H13 C12 1.1028 C11 109.2276 O11 -121.6416
            H14 C11 1.1025 C12 115.6994 H11  180.0000
            0 1
        '''

        atom_name = [
            'C', 'O', 'CB', 'HB1', 'HB2', 'HB3', 'HA'
        ]

        return textwrap.dedent(zmatrix), atom_name

