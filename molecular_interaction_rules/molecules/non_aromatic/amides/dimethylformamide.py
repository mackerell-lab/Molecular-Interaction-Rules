#!/usr/bin/env python3
#
# Lennard-Jones-Drill-2: Dimethylformamide
# ----------------------------------------

# Imports
# -------
import textwrap

class Dimethylformamide(object):

    def __init__(self):

        self.resi_name = 'DMF'

    def get_monomer_a_species(self):

        '''

        Get the Monomer A Species

        '''

        monomer_a_species = {
            'H1': self.get_monomer_a_amide_hydrogen_zmatrix(),
            'N1': self.get_monomer_a_nitrogen_zmatrix(),
        }

        return monomer_a_species

    def get_monomer_a_amide_hydrogen_zmatrix(self):

        zmatrix = '''\
          H11
          C12   H11   1.1026
          N11   C12   1.3795   H11  126.0143
          C11   N11   1.4499   C12  120.2747    H11  178.1853
          O11   C12   1.2241   N11  126.0143    C11   -2.4121
          C13   N11   1.4484   C12  121.2997    H11    6.6514
          H12   C13   1.0936   N11  109.0238    C11   57.4875
          H13   C13   1.0935   N11  111.1224    C11  177.9018
          H14   C13   1.0938   N11  108.9572    C11  -62.0719
          H15   C11   1.0938   N11  108.9572    C12  132.5108
          H16   C11   1.0939   N11  108.9572    C12 -108.1415
          H17   C11   1.0932   N11  111.1224    C12   12.0764
          X11   H11   1.0000   C12   90.0000    N11  180.0000
          0 1
        '''

        atom_name = [
          'HA', 'C', 'N', 'CC', 'O', 'CT', 'HC1', 'HC2', 'HC3', 'HT1', 'HT2', 'HT3'
        ]

        return zmatrix, atom_name

    def get_monomer_a_nitrogen_zmatrix(self):

        zmatrix = '''\
          N11
          C12   N11   1.3795
          H11   C12   1.1026   N11   126.0143
          C11   N11   1.4499   C12  120.2747    H11  178.1853
          O11   C12   1.2241   N11  126.0143    C11   -2.4121
          C13   N11   1.4484   C12  121.2997    H11    6.6514
          H12   C13   1.0936   N11  109.0238    C11   57.4875
          H13   C13   1.0935   N11  111.1224    C11  177.9018
          H14   C13   1.0938   N11  108.9572    C11  -62.0719
          H15   C11   1.0938   N11  108.9572    C12  132.5108
          H16   C11   1.0939   N11  108.9572    C12 -108.1415
          H17   C11   1.0932   N11  111.1224    C12   12.0764
          X11   N11   1.0000   C12   90.0000    H11  180.0000
          0 1
        '''

        atom_name = [
          'N', 'C', 'HA', 'CC', 'O', 'CT', 'HC1', 'HC2', 'HC3', 'HT1', 'HT2', 'HT3'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_b_species(self):

        '''

        Get the Monomer B Species

        '''

        monomer_b_species = {
            'H1': self.get_monomer_b_hydrogen_zmatrix(),
            # 'N1': self.get_monomer_b_nitrogen_zmatrix()
        }

        return monomer_b_species

    def get_monomer_b_nitrogen_zmatrix(self):

        zmatrix = '''\
          N21   :1  DISTANCE   :2  ANGLE        :3   90.00000
          C22   N21   1.3795   :1  180.0000     :2   DIHEDRAL
          H21   C22   1.1026   N21   126.0143   :1     0.0000
          C21   N21   1.4499   C22  120.2747    H21  178.1853
          O21   C22   1.2241   N21  126.0143    C21   -2.4121
          C23   N21   1.4484   C22  121.2997    H21    6.6514
          H22   C23   1.0936   N21  109.0238    C21   57.4875
          H23   C23   1.0935   N21  111.1224    C21  177.9018
          H24   C23   1.0938   N21  108.9572    C21  -62.0719
          H25   C21   1.0938   N21  108.9572    C22  132.5108
          H26   C21   1.0939   N21  108.9572    C22 -108.1415
          H27   C21   1.0932   N21  111.1224    C22   12.0764
          X21   N21   1.0000   C22   90.0000    H21  180.0000
          0 1
        '''

        atom_name = [
          'N', 'C', 'HA', 'CC', 'O', 'CT', 'HC1', 'HC2', 'HC3', 'HT1', 'HT2', 'HT3'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_b_hydrogen_zmatrix(self):

        zmatrix = '''\
            H21    :1   DISTANCE   :2  ANGLE       :3   DIHEDRAL
            X21  H21    1.0000     :1  90.0000     :2    0.0000
            C22   H21   1.1026   X21  90.0000      :2  180.0000
            N21   C22   1.3795   H21  126.0143     :1  180.0000
            C21   N21   1.4499   C22  120.2747    H21  178.1853
            O21   C22   1.2241   N21  126.0143    C21   -2.4121
            C23   N21   1.4484   C22  121.2997    H21    6.6514
            H22   C23   1.0936   N21  109.0238    C21   57.4875
            H23   C23   1.0935   N21  111.1224    C21  177.9018
            H24   C23   1.0938   N21  108.9572    C21  -62.0719
            H25   C21   1.0938   N21  108.9572    C22  132.5108
            H26   C21   1.0939   N21  108.9572    C22 -108.1415
            H27   C21   1.0932   N21  111.1224    C22   12.0764
            0 1
        '''

        atom_name = [
          'HA', 'C', 'N', 'CC', 'O', 'CT', 'HC1', 'HC2', 'HC3', 'HT1', 'HT2', 'HT3'
        ]

        return textwrap.dedent(zmatrix), atom_name

