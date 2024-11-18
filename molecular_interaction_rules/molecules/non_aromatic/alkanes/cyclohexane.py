# Imports
# -------

import textwrap

class Cyclohexane(object):

    def __init__(self):

        self.resi_name = 'CHEX'

    def get_monomer_a_species(self):

        '''

        Get the Monomer A Species

        '''

        monomer_a_species = {
            'H1': self.get_chair_hydrogen_conformation()
        }

        return monomer_a_species

    def get_monomer_b_species(self):

        monomer_b_species = {
            'H1': self.get_monomer_b_chair_hydrogen_conformation(),
        }

        return monomer_b_species

    def get_chair_carbon_conformation(self):

        zmatrix = '''\
            C11
            C12 C11 1.5369
            C13 C12 1.5369 C11 110.9384
            C14 C13 1.5369 C12 110.9384 C11  -56.214
            C15 C13 2.5322 C12  90.0000 C11  -28.107
            C16 C13 2.9621 C12  58.7453 C11  -28.1071
            H11 C14 1.1035 C13 110.1108 C12 178.8057
            H12 C14 1.1064 C13 110.1108 C12  -64.0225
            H13 C15 1.1064 C14 109.1109 C13   64.0225
            H14 C15 1.1064 C14 109.1109 C13 -178.8056
            H15 C16 1.1064 C11 109.1108 C12   64.0224
            H16 C16 1.1064 C11 110.3361 C12  -178.8058
            H17 C13 1.1035 C12 110.3361 C11  -178.8057
            H18 C13 1.1064 C12 109.1108 C11    64.0225
            H19 C12 1.1064 C13 109.1108 C14    64.0224
            H110 C12 1.1064 C13 110.3362 C14 -178.8058
            H111 C11 1.1064 C12 110.3362 C13  -64.0224
            H112 C11 1.1064 C12 110.3361 C13  178.8058
            X11  C11 1.0000 C12  90.0000 C13  180.0000
            0 1
        '''

        atom_name = [
          'C1', 'C2', 'C3', 'C4', 'C5', 'C6',
          'H4A', 'H4B', 'H5A', 'H5B', 'H6A', 'H6B',
          'H1A', 'H1B', 'H2A', 'H2B', 'H3A', 'H3B'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_chair_hydrogen_conformation(self):

        zmatrix = '''\
            H11
            C11 H11 1.1064
            C12 C11 1.5369 H11 110.3361
            C13 C12 1.5369 C11 110.9384 H11 178.8058
            C14 C13 1.5369 C12 110.9384 C11  -56.214
            C15 C13 2.5322 C12  90.0000 C11  -28.107
            C16 C13 2.9621 C12  58.7453 C11  -28.1071
            H11 C14 1.1035 C13 110.1108 C12 178.8057
            H12 C14 1.1064 C13 110.1108 C12  -64.0225
            H13 C15 1.1064 C14 109.1109 C13   64.0225
            H14 C15 1.1064 C14 109.1109 C13 -178.8056
            H15 C16 1.1064 C11 109.1108 C12   64.0224
            H16 C16 1.1064 C11 110.3361 C12  -178.8058
            H17 C13 1.1035 C12 110.3361 C11  -178.8057
            H18 C13 1.1064 C12 109.1108 C11    64.0225
            H19 C12 1.1064 C13 109.1108 C14    64.0224
            H110 C12 1.1064 C13 110.3361 C14 -178.8058
            H111 C11 1.1064 C12 110.3361 C13  -64.0224
            X11  H11 1.0000 C11  90.0000 C12  180.0000
            0 1
        '''

        atom_name = [
          'H3B', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6',
          'H4A', 'H4B', 'H5A', 'H5B', 'H6A', 'H6B',
          'H1A', 'H1B', 'H2A', 'H2B', 'H3A',
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_b_chair_hydrogen_conformation(self):

        zmatrix = '''\
            H21  :1 DISTANCE  :2 ANGLE :3  DIHEDRAL
            X21 H21 1.0000  :1  90.0000   :2   0.0000
            C21 H21 1.1064 X21  90.0000   :1  180.0000
            C22 C21 1.5369 H21 110.3361   :1  180.0000
            C23 C22 1.5369 C21 110.9384 H21 178.8058
            C24 C23 1.5369 C22 110.9384 C21  -56.214
            C25 C23 2.5322 C22  90.0000 C21  -28.107
            C26 C23 2.9621 C22  58.7453 C21  -28.1071
            H21 C24 1.1035 C23 110.1108 C22 178.8057
            H22 C24 1.1064 C23 110.1108 C22  -64.0225
            H23 C25 1.1064 C24 109.1109 C23   64.0225
            H24 C25 1.1064 C24 109.1109 C23 -178.8056
            H25 C26 1.1064 C21 109.1108 C22   64.0224
            H26 C26 1.1064 C21 110.3361 C22  -178.8058
            H27 C23 1.1035 C22 110.3361 C21  -178.8057
            H28 C23 1.1064 C22 109.1108 C21    64.0225
            H29 C22 1.1064 C23 109.1108 C24    64.0224
            H210 C22 1.1064 C23 110.3361 C24 -178.8058
            H211 C21 1.1064 C22 110.3361 C23  -64.0224
            0 1
          '''

        atom_name = [
          'H3B', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6',
          'H4A', 'H4B', 'H5A', 'H5B', 'H6A', 'H6B',
          'H1A', 'H1B', 'H2A', 'H2B', 'H3A',
        ]

        return textwrap.dedent(zmatrix), atom_name

