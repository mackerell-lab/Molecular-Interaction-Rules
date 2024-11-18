#!/usr/bin/env python3
#
# Lennard-Jones-Drill-2: Dimethyl Phosphate
# -----------------------------------------

# Imports
# -------

import textwrap
#!/usr/bin/env python3
#
# Lennard-Jones-Drill-2: Dimethyl Phosphate
# -----------------------------------------

# Imports
# -------

import textwrap

class DimethylPhosphate(object):

    def __init__(self):

        self.resi_name = 'DMEP'

    def get_monomer_a_species(self):

        '''

        Get the Monomer A Species

        '''

        monomer_a_species = {
            'P1': self.get_double_bond_oxygen_lone_pair_donor()
        }

        return monomer_a_species

    def get_double_bond_oxygen_lone_pair_donor(self):

        zmatrix = '''\
            P11
            O11 P11 1.6437
            O12 P11 1.6437 O11 97.8486
            H11 O12 0.9705 P11 108.6226 O11 146.0940
            O13 P11 1.6286 O12 104.4218 H11 -106.7457
            C11 O13 1.4523 P11 117.8616 O11 51.1114
            H12 C11 1.0992 O13 110.2301 P11 -60.9713
            H13 C11 1.0992 O13 110.2301 P11 60.9714
            H14 C11 1.0957 O13 105.5916 P11 -179.9999
            O14 P11 1.5016 O12 117.4240 H11 19.5590
            H15 O11 0.9705 P11 108.6223 O12 -146.0925
            X11 O14 1.0000 P11 90.0000 O13 180.0000
            0 1
        '''

        atom_name = [
            'O3', 'P1', 'O4', ''
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_b_species(self):

      '''

      Get the Monomer B Species

      '''

      monomer_b_species = {
          'P1': self.get_monomer_b_oxygen_zmatrix()
      }

      return monomer_b_species

    def get_monomer_b_oxygen_zmatrix(self):

        zmatrix = '''\
            O21   :1  DISTANCE  :2  ANGLE  :3   DIHEDRAL
            P21  O21  1.5085   :1   180.0000  :2  180.0000
            O22  P21  1.6377  O21  115.2694   :1    0.0000
            H21  O22  0.9709  P21  108.4278  O21  -21.0969
            O23  P21  1.6251  O22  101.1193  H21  107.3959
            O24  P21  1.6311  O22  102.5485  H21 -148.9298
            C21  O23  1.4551  P21  115.2049  O21  -50.3827
            C22  O24  1.4542  P21  115.4507  O21  -44.3740
            C23  O24  1.4542  P21  115.4507  O21  -44.3740
            H22  C21  1.0994  O23  109.8941  P21  -68.8607
            H23  C21  1.0994  O23  109.8941  P21   53.1472
            H24  C21  1.0994  O23  109.8941  P21  172.3540
            H25  C22  1.0994  O24  109.8941  P21  170.8612
            H26  C22  1.0994  O24  109.8941  P21   51.3387
            H27  C22  1.0994  O24  109.8941  P21  -70.4974
            X21  O21  1.0000  P21   90.0000  O22  180.0000
            0 1
          '''

        atom_name = [
          ''
        ]

        return textwrap.dedent(zmatrix), atom_name

