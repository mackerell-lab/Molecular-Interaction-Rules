#!/usr/bin/env python3
#
# Lennard-Jones-Drill-2: Prolineamide
# -----------------------------------

# Imports
# -------
import textwrap

class Prolineamide(object):

    def __init__(self):

        self.resi_name = 'PNH1'

    def get_monomer_a_species(self):

        '''

        Get the Monomer A Species

        '''

        monomer_a_species = {
          'H1': self.get_monomer_a_amide_nitrogen(),
          'H2': self.get_monomer_a_sp3_carbon_one_hydrogen(),
        }

        return monomer_a_species

    def get_monomer_a_amide_nitrogen(self):

        zmatrix = '''\
                H11
                N12      H11    1.0236
                C12      N12    1.4911     H11   107.8309
                C11      C12     1.5371     N12  109.2871     H11  153.7244
                N11      C11    1.3773     C12    113.2503    N12  -36.7063
                C13      N12    1.4808      C12  103.0620      C11  -92.1041
                H12      C13    1.1011      N12  110.2057      C12 -163.3153
                H13      C13    1.1036      N12  107.5679     C12   79.3090
                C14      C13    1.5435      N12  106.3270      H11   74.8116
                H14      C14    1.1029      C13  109.2528      N12  -84.9363
                H15      C14    1.1019      C13  113.1750      N12  155.3921
                C15      C12    1.5530      C11  111.8529      N11 -156.4619
                H16      C15    1.1021      C12  110.4009     C11  -10.1436
                H17      C15    1.1008      C12  111.0245      C11 -128.7878
                H18      C12    1.1020      C11  106.0926      N11   80.4512
                O11      C11    1.2355      C12  122.3919      N12  148.7052
                H19      N11    1.0137      C11  116.2255      C12  171.7158
                H20      N11    1.0180      C11  116.8170      C12   23.3230
                X11      H11     1.0000     N12  90.0000       C12   180.0000
                0 1
        '''

        atom_name = [
          'H9', 'N2', 'C5', 'C1', 'N1', 'C4', 'H8', 'C3', 'C2', 'O', 'H6', 'H7', 'H4', 'H5', 'H2', 'H3', 'H1', 'H10'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_a_sp3_carbon_one_hydrogen(self):

        zmatrix = '''\
            H11
            C12  H11  1.0846
            C11  C12  1.5300  H11  104.5232
            N11  C11  1.3445  C12  115.3149  H11  -97.9238
            N12  C12  1.4641  C11  112.2071  N11   20.4696
            C13  N12  1.4615  C12  105.8259  C11   99.2627
            H12  N12  1.0010  C12  109.4950  C11 -143.0365
            C14  C13  1.5322  N12  105.5222  H12  -81.1190
            C15  C12  1.5300  C11  111.9567  N11  141.5078
            O11  C11  1.2019  C12  121.1157  N12 -162.3860
            H13  C13  1.0846  N12  108.5988  C12  -81.7117
            H14  C13  1.0846  N12  110.5262  C12  160.4969
            H15  C14  1.0846  C13  112.8941  N12 -156.2638
            H16  C14  1.0846  C13  112.8941  N12   83.8645
            H17  C15  1.0846  C12  110.2896  C11   -0.3212
            H18  C15  1.0846  C12  111.1976  C11  118.1844
            H19  N11  0.9942  C11  120.5945  C12 -179.6150
            H20  N11  0.9942  C11  120.5945  C12   -8.4508
            X11  H11  1.0000  C12   90.0000  C11  180.0000
            0 1
        '''

        atom_name = [
          'H1', 'C1', 'C5', 'N2', 'N1', 'C4', 'H8', 'C3', 'C2', 'O', 'H6', 'H7', 'H4', 'H5', 'H2', 'H3', 'H9', 'H10'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_b_species(self):

        '''

        Get the Monomer B Species

        '''

        monomer_b_species = {
            'H1': self.get_monomer_b_nitrogen_hydrogen_zmatrix(),
            'H2': self.get_monomer_b_hydrogen_zmatrix()
        }

        return monomer_b_species

    def get_monomer_b_nitrogen_hydrogen_zmatrix(self):

        zmatrix = '''\
            H21       :1   DISTANCE      :2  ANGLE          :3   DIHEDRAL
            X21      H21    1.0000       :1  90.0000        :2    0.0000
            N22      H21    1.0236      X21   90.0000       :2  180.0000
            C22      N22    1.4911      H21   107.8309      :1  180.0000
            C21      C22     1.5371     N22  109.2871      H21  153.7244
            N21      C21    1.3773      C22  113.2503      N22  -36.7063
            C23      N22    1.4808      C22  103.0620      C21  -92.1041
            H22      C23    1.1011      N22  110.2057      C22 -163.3153
            H23      C23    1.1036      N22  107.5679      C22   79.3090
            C24      C23    1.5435      N22  106.3270      H21   74.8116
            H24      C24    1.1029      C23  109.2528      N22  -84.9363
            H25      C24    1.1019      C23  113.1750      N22  155.3921
            C25      C22    1.5530      C21  111.8529      N21 -156.4619
            H26      C25    1.1021      C22  110.4009      C21  -10.1436
            H27      C25    1.1008      C22  111.0245      C21 -128.7878
            H28      C22    1.1020      C21  106.0926      N21   80.4512
            O21      C21    1.2355      C22  122.3919      N22  148.7052
            H29      N21    1.0137      C21  116.2255      C22  171.7158
            H30      N21    1.0180      C21  116.8170      C22   23.3230
         0 1
        '''

        atom_name = [
          'H9', 'N2', 'C5', 'C1', 'N1', 'C4', 'H8', 'C3', 'C2', 'O', 'H6', 'H7', 'H4', 'H5', 'H2', 'H3', 'H1', 'H10'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_b_hydrogen_zmatrix(self):

        zmatrix = '''\
            H21    :1   DISTANCE   :2  ANGLE   :3   DIHEDRAL
            X21  H21    1.0000     :1  90.0000   :2    0.0000
            C22  H21  1.0846  X21  90.0000       :2  180.0000
            C21  C22  1.5300  H21  104.5232      :1  180.0000
            N21  C21  1.3445  C22  115.3149  H21  -97.9238
            N22  C22  1.4641  C21  112.2071  N21   20.4696
            C23  N22  1.4615  C22  105.8259  C21   99.2627
            H22  N22  1.0010  C22  109.4950  C21 -143.0365
            C24  C23  1.5322  N22  105.5222  H22  -81.1190
            C25  C22  1.5300  C21  111.9567  N21  141.5078
            O21  C21  1.2019  C22  121.1157  N22 -162.3860
            H23  C23  1.0846  N22  108.5988  C22  -81.7117
            H24  C23  1.0846  N22  110.5262  C22  160.4969
            H25  C24  1.0846  C23  112.8941  N22 -156.2638
            H26  C24  1.0846  C23  112.8941  N22   83.8645
            H27  C25  1.0846  C22  110.2896  C21   -0.3212
            H28  C25  1.0846  C22  111.1976  C21  118.1844
            H29  N21  0.9942  C21  120.5945  C22 -179.6150
            H30  N21  0.9942  C21  120.5945  C22   -8.4508
            0 1
        '''

        atom_name = [
          'H1', 'C1', 'C5', 'N2', 'N1', 'C4', 'H8', 'C3', 'C2', 'O', 'H6', 'H7', 'H4', 'H5', 'H2', 'H3', 'H9', 'H10'
        ]

        return textwrap.dedent(zmatrix), atom_name

