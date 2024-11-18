#!/usr/bin/env python3
#
# Lennard-Jones-Drill-2: DimethylTrithiocarbonate
# -----------------------------------------------

# Imports
# -------

import textwrap

class DimethylTrithiocarbonate(object):

    def __init__(self):

        self.resi_name = 'DMTT'

    def get_monomer_a_species(self):

        '''

        Get the Monomer A Species

        '''

        monomer_a_species = {
            'S1': self.get_conjoining_sulphur(),
            'S2': self.get_carbonyl_sulphur(),
            'C1': self.get_carbon_sulphur_zmatrix()
        }

        return monomer_a_species

    def get_conjoining_sulphur(self):

        zmatrix = '''\
          S11
          C11  S11 1.7783
          S12  C11 1.7783  S11  118.9556
          C12  S12 1.8317  C11  102.3499  S11   42.4811
          H11  C12 1.0964  S12  113.1207  C11  -70.8872
          H12  C12 1.1006  S12  105.1031  C11  169.7883
          H13  C12 1.0981  S12  109.2015  C11   52.4264
          S13  C11 1.6514  S12  120.5202  C12 -137.5173
          C13  S11 1.8317  C11  102.3421  S12   42.5371
          H14  C13 1.0981  S11  109.2007  C11   52.4852
          H15  C13 1.0981  S11  105.1053  C11  169.8482
          H16  C13 1.0981  S11  113.1186  C11  -70.8254
          X11  S11 1.0000 C11  135.0000  C12  -90.0000
          0 1
        '''

        atom_name = [
          'S1', 'C', 'S2', 'C2', 'H21', 'H22', 'H23', 'S', 'C1', 'H11', 'H12', 'H13'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_carbon_sulphur_zmatrix(self):

        zmatrix = '''\
          C11
          S11  C11 1.7783
          S12  C11 1.7783  S11  118.9556
          C12  S12 1.8317  C11  102.3499  S11   42.4811
          H11  C12 1.0964  S12  113.1207  C11  -70.8872
          H12  C12 1.1006  S12  105.1031  C11  169.7883
          H13  C12 1.0981  S12  109.2015  C11   52.4264
          S13  C11 1.6514  S12  120.5202  C12 -137.5173
          C13  S11 1.8317  C11  102.3421  S12   42.5371
          H14  C13 1.0964  S11  109.2007  C11   52.4852
          H15  C13 1.0964  S11  105.1053  C11  169.8482
          H16  C13 1.0964  S11  113.1186  C11  -70.8254
          X11  C11 1.0000  S13   90.0000  S11   90.0000
          0 1
        '''

        atom_name = [
          'C','S1','S2', 'C2', 'H21', 'H22', 'H23', 'S', 'C1', 'H11', 'H12', 'H13'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_carbonyl_sulphur(self):

        zmatrix = '''\
          S13
          C11  S13 1.6514
          S12  C11 1.7783  S13  120.5202
          C12  S12 1.8317  C11  102.3499  S13 -137.5173
          S11  C11 1.7783  S12  118.9556  C12   45.0000
          C13  S11 1.8317  C11  102.3421  S12   42.5371
          H11  C12 1.0964  S12  113.1207  C11  -70.8872
          H12  C12 1.1006  S12  105.1031  C11  169.7883
          H13  C12 1.0981  S12  109.2015  C11   52.4264
          H14  C13 1.0981  S11  109.2007  C11   52.4852
          H15  C13 1.1006  S11  105.1053  C11  169.8482
          H16  C13 1.0964  S11  113.1186  C11  -70.8254
          X11  S13 1.0000  C11   90.0000  S12   90.0000
          0 1
        '''

        atom_name = [
          'S', 'C', 'S1', 'C2', 'S2', 'C1', 'H21', 'H22', 'H23', 'H11', 'H12', 'H13'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_b_species(self):

        '''

        Get the Monomer B Species

        '''

        monomer_b_species = {
            'S1': self.get_monomer_b_conjoining_sulphur(),
            'S2': self.get_monomer_b_carbonyl_sulphur(),
        }

        return monomer_b_species

    def get_monomer_b_conjoining_sulphur(self):

        zmatrix = '''\
            S21   :1  DISTANCE   :2  ANGLE    :3    DIHEDRAL
            C21  S21 1.7783      :1  180.0000    :2   180.0000
            S22  C21 1.7783  S21  118.9556       :1     0.0000
            C22  S22 1.8317  C21  102.3499  S21   42.4811
            H21  C22 1.0964  S22  113.1207  C21  -70.8872
            H22  C22 1.1006  S22  105.1031  C21  169.7883
            H23  C22 1.0981  S22  109.2015  C21   52.4264
            S23  C21 1.6514  S22  120.5202  C22 -137.5173
            C23  S21 1.8317  C21  102.3421  S22   42.5371
            H24  C23 1.0981  S21  109.2007  C21   52.4852
            H25  C23 1.1006  S21  105.1053  C21  169.8482
            H26  C23 1.0964  S21  113.1186  C21  -70.8254
            X21  S21 1.0000  C21  135.0000  C22  -90.0000
            0 1
          '''

        atom_name = [
          'S1', 'C', 'S2', 'C2', 'H21', 'H22', 'H23', 'S', 'C1', 'H11', 'H12', 'H13'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_b_carbonyl_sulphur(self):

        zmatrix = '''\
            S23   :1  DISTANCE   :2  ANGLE    :3    DIHEDRAL
            C21  S23 1.6514      :1  180.0000    :2   180.0000
            S22  C21 1.7783  S23  120.5202       :1     0.0000
            C22  S22 1.8317  C21  102.3499  S23 -137.5173
            S21  C21 1.7783  S22  118.9556  C22   45.0000
            C23  S21 1.8317  C21  102.3421  S22   42.5371
            H21  C22 1.0964  S22  113.1207  C21  -70.8872
            H22  C22 1.1006  S22  105.1031  C21  169.7883
            H23  C22 1.0981  S22  109.2015  C21   52.4264
            H24  C23 1.0981  S21  109.2007  C21   52.4852
            H25  C23 1.1006  S21  105.1053  C21  169.8482
            H26  C23 1.0964  S21  113.1186  C21  -70.8254
            X21  S23 1.0000  C21   90.0000  S22   90.0000
            0 1
          '''

        atom_name = [
          'S', 'C', 'S1', 'C2', 'S2', 'C1', 'H21', 'H22', 'H23', 'H11', 'H12', 'H13'
        ]

        return textwrap.dedent(zmatrix), atom_name

