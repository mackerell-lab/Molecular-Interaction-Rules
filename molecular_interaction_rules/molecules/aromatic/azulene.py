#!/usr/bin/env python3
#
# Lennard-Jones-Drill-2: Azulene
# ------------------------------

# Imports
# -------
import textwrap

class Azulene(object):

    def __init__(self):

        self.resi_name = 'azul'

    def get_monomer_a_species(self):

        '''

        Get the Monomer A Species

        '''

        monomer_a_species = {
            'H1': self.monomer_a_hydrogen_zmatrix(),
            'RC1': self.monomer_a_aromatic_zmatrix()
        }

        return monomer_a_species

    def get_monomer_b_species(self):

        monomer_b_species = {
            'H1': self.monomer_b_shaped_hydrogen_zmatrix(),
            'RC1': self.monomer_b_pi_stack_zmatrix()
        }

        return monomer_b_species

    def monomer_a_hydrogen_zmatrix(self):

        zmatrix = '''\
            H11
            C11   H11  1.0936
            C12   C11  1.4090   H11  115.5882
            C13   C12  1.4090   C11  128.9389    H11  180.0000
            C14   C13  1.4090   C12  128.9389    C11   -0.0000
            C15   C14  1.4090   C13  125.6525    C12 -180.0000
            C16   C15  1.4090   C14  108.3817    C13 -180.0000
            H12   C15  1.0921   C14  125.6525    C13    0.0000
            C17   C16  1.4090   C15  109.9255    H12  -180.0000
            C18   C13  2.6128   C12  101.3243    C11   -0.0000
            C19   C18  1.4090   C14  127.6958    C13   -0.0000
            C20   C11  1.4090   C12  129.6993    C13   -0.0000
            H13   C16  1.0936   C15  125.6525    C14 -180.0000
            H14   C17  1.0936   C16  126.5100    C15 -180.0000
            H15   C19  1.0936   C18  115.5882    C14 -180.0000
            H16   C20  1.0936   C11  115.5882    C12 -180.0000
            H17   C13  1.0936   C12  116.0586    C11 -180.0000
            H18   C12  1.0936   C13  115.5882    C14 -180.0000
            X11   H11  1.0000   C11   90.0000    C12  180.0000
            0 1
        '''

        atom_name = [
          'H5', 'C7', 'C6', 'C5', 'C4', 'C3', 'C2', 'H3', 'C1', 'C10', 'C9', 'C8', 'H2', 'H1', 'H9', 'H8', 'H7', 'H6'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def monomer_a_aromatic_zmatrix(self):

        zmatrix = '''\
          X11
          C11   X11  1.4000
          C12   C11  1.4090   X11   64.2850
          C13   C12  1.4090   C11  128.9389    X11    0.0000
          C14   C13  1.4090   C12  128.9389    C11   -0.0000
          C15   C14  1.4090   C13  125.6525    C12 -180.0000
          C16   C15  1.4090   C14  108.3817    C13 -180.0000
          H11   C15  1.0921   C14  125.6525    C13    0.0000
          C17   C16  1.4090   C15  109.9255    H11 -180.0000
          C18   C13  2.6128   C12  101.3243    C11   -0.0000
          C19   C18  1.4090   C14  127.6958    C13   -0.0000
          C20   C11  1.4090   C12  129.6993    C13   -0.0000
          H12   C16  1.0936   C15  125.6525    C14 -180.0000
          H13   C17  1.0936   C16  126.5100    C15 -180.0000
          H14   C19  1.0936   C18  115.5882    C14 -180.0000
          H15   C20  1.0936   C11  115.5882    C12 -180.0000
          H16   C13  1.0936   C12  116.0586    C11 -180.0000
          H17   C12  1.0936   C13  115.5882    C14 -180.0000
          H18   C11  1.0936   C12  115.5882    C13 -180.0000
          0 1
        '''

        atom_name = [
          'C7', 'C6', 'C5', 'C4', 'C3', 'C2', 'H3', 'C1', 'C10', 'C9', 'C8', 'H2', 'H1', 'H9', 'H8', 'H7', 'H6', 'H5'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def monomer_b_shaped_hydrogen_zmatrix(self):

        '''
        For replusion the angle is 180, and the dihedral is 180
        '''

        zmatrix = '''\
                H21    :1  DISTANCE  :2   ANGLE    :3   90.0000
                X21   H21  1.0000    :1   90.0000     :2   DIHEDRAL
                C21   H21  1.0756   X21   90.0000     :1  180.0000
                C22   C21  1.3774   H21  120.0000    X21    0.0000
                C23   C22  1.4090   C21  128.9389    H21  180.0000
                C24   C23  1.4090   C22  128.9389    C21   -0.0000
                C25   C24  1.4090   C23  125.6525    C22 -180.0000
                C26   C25  1.4090   C24  108.3817    C23 -180.0000
                H22   C25  1.0921   C24  125.6525    C23    0.0000
                C27   C26  1.4090   C25  109.9255    H22  -180.0000
                C28   C23  2.6128   C22  101.3243    C21   -0.0000
                C29   C28  1.4090   C24  127.6958    C23   -0.0000
                C30   C21  1.4090   C22  129.6993    C23   -0.0000
                H23   C26  1.0936   C25  125.6525    C24 -180.0000
                H24   C27  1.0936   C26  126.5100    C25 -180.0000
                H25   C29  1.0936   C28  115.5882    C24 -180.0000
                H26   C30  1.0936   C21  115.5882    C22 -180.0000
                H27   C23  1.0936   C22  116.0586    C21 -180.0000
                H28   C22  1.0936   C23  115.5882    C24 -180.0000
                0 1
            '''

        atom_name = [
          'H5', 'C7', 'C6', 'C5', 'C4', 'C3', 'C2', 'H3', 'C1', 'C10', 'C9', 'C8', 'H2', 'H1', 'H9', 'H8', 'H7', 'H6'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def monomer_b_pi_stack_zmatrix(self):

        zmatrix = '''\
            X21    :1  DISTANCE  :2   ANGLE     :3   90.0000
            C21   X21  1.3940    :1   90.0000     :2  180.0000
            C22   C21  1.3774   X21   60.0000     :1   90.0000
            C23   C22  1.3774   C21  -120.0000   X21   DIHEDRAL
            C24   C23  1.4090   C22  128.9389    C21   -0.0000
            C25   C24  1.4090   C23  125.6525    C22 -180.0000
            C26   C25  1.4090   C24  108.3817    C23 -180.0000
            H21   C25  1.0921   C24  125.6525    C23    0.0000
            C27   C26  1.4090   C25  109.9255    H21 -180.0000
            C28   C23  2.6128   C22  101.3243    C21   -0.0000
            C29   C28  1.4090   C24  127.6958    C23   -0.0000
            C30   C21  1.4090   C22  129.6993    C23   -0.0000
            H22   C26  1.0936   C25  125.6525    C24 -180.0000
            H23   C27  1.0936   C26  126.5100    C25 -180.0000
            H24   C29  1.0936   C28  115.5882    C24 -180.0000
            H25   C30  1.0936   C21  115.5882    C22 -180.0000
            H26   C23  1.0936   C22  116.0586    C21 -180.0000
            H27   C22  1.0936   C23  115.5882    C24 -180.0000
            H28   C21  1.0936   C22  115.5882    C23 -180.0000
            0 1
        '''

        atom_name = [
          'C7', 'C6', 'C5', 'C4', 'C3', 'C2', 'H3', 'C1', 'C10', 'C9', 'C8', 'H2', 'H1', 'H9', 'H8', 'H7', 'H6', 'H5'
        ],

        return textwrap.dedent(zmatrix), atom_name

