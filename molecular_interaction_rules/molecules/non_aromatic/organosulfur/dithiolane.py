#!/usr/bin/env python3
#
# Lennard-Jones-Drill-2: Dithiolane
# ---------------------------------

# Imports
# -------

import textwrap

class Dithiolane(object):

    def __init__(self):

        self.resi_name = 'MESH'

    def get_monomer_a_species(self):

        '''

        Get the Monomer A Species

        '''

        monomer_a_species = {
        }

        return monomer_a_species

    def get_monomer_a(self):

        zmatrix = '''\
        '''

        atom_name = [
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_b_species(self):

        monomer_b_species = {
        }

        return monomer_b_species

    def get_monomer_b(self):

        zmatrix = '''\
        '''

        atom_name = [
        ]

        return textwrap.dedent(zmatrix), atom_name
