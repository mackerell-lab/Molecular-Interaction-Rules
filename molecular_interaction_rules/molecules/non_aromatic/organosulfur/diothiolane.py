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

    def get_sulphur_hetereoatom(self):

        zmatrix = '''\
            S11
            S12  S11  2.1492
            C11  S12  1.8415    S11   94.9681
            C12  C11  1.5287    S12  108.1008    S11  -28.2456
            H11  C12  1.1043    C11  110.1289    S12  173.1128
            H12  C12  1.1009    C11  109.8789    S12  -66.6192
            C13  S11  1.8415    S12   94.9680    C11   -0.0003
            H13  C13  1.1011    S11  106.2198    S12  148.9234
            H14  C13  1.1009    S11  109.4742    S12  -92.8470
            H15  C11  1.1011    S12  106.2199    S11 -148.9227
            H16  C11  1.1009    S12  109.4741    S11   92.8477
        '''

        atom_name = []
        return textwrap.dedent(zmatrix), atom_name
      
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
