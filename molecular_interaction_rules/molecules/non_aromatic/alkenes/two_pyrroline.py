# Imports
# -------

import textwrap

class TwoPyrroline(object):

    def __init__(self):

        self.resi_name = '2PRL'

    def get_monomer_a_species(self):

        '''

        Get the Monomer A Z-Matrices

        '''

        monomer_a_species = {
            'H1': self.get_monomer_a_nitrogen_zmatrix()
        }

        return monomer_a_species

    def get_monomer_a_nitrogen_zmatrix(self):

        zmatrix = '''\
            H11
            N11  H11  0.9989
            C11  N11  1.4049  H11  114.3393
            C12  C11  1.3196  N11  113.4180  H11  142.1118
            C13  C12  1.5133  C11  109.0047  N11   -0.2544
            H12  C13  1.0881  C12  111.0556  C11  103.9786
            H13  C13  1.0851  C12  113.9511  C11 -134.8276
            C14  N11  1.4675  C11  105.9046  C12   15.8159
            H14  C14  1.0886  N11  110.6776  C11   94.1438
            H15  C14  1.0824  N11  110.8943  C11 -145.7516
            H16  C12  1.0721  C11  125.9109  N11  174.3339
            H17  C11  1.0737  C12  127.6615  C13 -176.6513
            X11  H11  1.0000  N11   90.0000  C11    0.0000
            0 1
        '''

        atom_name = [
            'H1', 'N1', 'C2', 'C3', 'C4', 'H41', 'H42', 'C5', 'H51', 'H52', 'H3', 'H2'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_b_species(self):

        '''

        Get the Monomer A Z-Matrices

        '''

        monomer_b_species = {
        }

        return monomer_b_species

