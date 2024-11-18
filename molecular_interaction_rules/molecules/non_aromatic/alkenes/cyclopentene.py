# Imports
# -------

import textwrap

class Cyclopentene(object):
    
    def __init__(self):

        self.resi_name = 'CYPE'

    def get_monomer_a_species(self):

        '''

        Get the Monomer A Z-Matrices

        '''

        monomer_a_species = {
          'H1': self.get_monomer_a_sp2_zmatrix(),
          'H2': self.get_monomer_a_sp3_zmatrix()
        }

        return monomer_a_species

    def get_monomer_b_species(self):

      '''

      Get the Monomer A Z-Matrices

      '''

      monomer_b_species = {
        'H1': self.get_monomer_b_sp2_zmatrix(),
        'H2': self.get_monomer_b_sp3_zmatrix()
      }

      return monomer_b_species

    def get_monomer_a_sp2_zmatrix(self):

      zmatrix = '''\
          H11
          C13   H11  1.0935
          C12   C13  1.5178   H11  123.7585
          C11   C12  1.5508   C13  102.4718   H11 -166.6366
          C14   C13  1.3545   C12  111.3707   C11   16.4054
          C15   C11  1.5508   C12  105.0372   C13  -25.6092
          H12   C14  1.1007   C13  124.7945   C12  176.9201
          H13   C15  1.1007   C11  111.3707   C12  -92.1016
          H14   C15  1.1007   C11  112.5420   C12  147.5060
          H15   C12  1.1007   C13  113.1858   C14  137.8623
          H16   C12  1.1007   C13  110.0838   C14 -102.2965
          H17   C11  1.1007   C12  108.7004   C13   90.5856
          H18   C11  1.1007   C12  113.1858   C13 -149.3354
          X11   H11  1.0000   C13   90.0000   C12  180.0000
          0 1
        '''

      atom_name = [
        'H2', 'C2', 'C3', 'C4', 'C1', 'C5', 'H1', 'H51', 'H52', 'H2', 'H31', 'H32', 'H42', 'H41',
      ]

      return textwrap.dedent(zmatrix), atom_name

    def get_monomer_b_sp2_zmatrix(self):

      zmatrix = '''\
          H21     :1    DISTANCE     :2  ANGLE        :3  DIHEDRAL
          X21     H21   1.0000       :1  90.0000        :2    0.0000
          C23    H21    1.0935      X21  90.0000        :2  180.0000
          C22    C23    1.5353      H21  117.8400       :1  180.0000
          C21   C22  1.5505   C23  102.4718   H21 -166.6366
          C24   C23  1.35468  C22  111.3707   C21   16.4054
          C25   C21  1.5505   C22  105.0372   C23  -25.6092
          H22   C24  1.1007   C23  124.7945   C22  176.9201
          H23   C25  1.1007   C21  111.3707   C22  -92.1016
          H24   C25  1.1007   C21  112.5420   C22  147.5060
          H25   C22  1.1007   C23  113.1858   C24  137.8623
          H26   C22  1.1007   C23  110.0838   C24 -102.2965
          H27   C21  1.1007   C22  108.7004   C23   90.5856
          H28   C21  1.1007   C22  113.1858   C23 -149.3354
          X21   H21  1.0000   C23   90.0000   C22  180.0000
          0 1
        '''

      atom_name = [
        'H2', 'C2', 'C3', 'C4', 'C1', 'C5', 'H1', 'H51', 'H52', 'H2', 'H31', 'H32', 'H42', 'H41',
      ]

      return textwrap.dedent(zmatrix), atom_name

    def get_monomer_a_sp3_zmatrix(self):

        zmatrix = '''\
          H11
          C11   H11  1.1007
          C12   C11  1.5505   H11  113.1858
          C13   C12  1.5178   C11  102.4718   H11 -149.3354
          C14   C13  1.3546   C12  111.3707   C11   16.4054
          C15   C11  1.5505   C12  105.0372   C13  -25.6092
          H12   C14  1.1007   C13  124.7945   C12  176.9201
          H13   C15  1.1007   C11  111.3707   C12  -92.1016
          H14   C15  1.1007   C11  112.5420   C12  147.5060
          H15   C13  1.1007   C12  123.7585   C11 -166.6366
          H16   C12  1.1007   C13  113.1858   C14  137.8623
          H17   C12  1.1007   C13  110.0838   C14 -102.2965
          H18   C11  1.1007   C12  108.7004   C13   90.5856
          X11   H11  1.0000   C11   90.0000   C12  180.0000
          0 1
        '''

        atom_name = [
          'H41', 'C4', 'C3', 'C2', 'C1', 'C5', 'H1', 'H51', 'H52', 'H2', 'H31', 'H32', 'H42'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def get_monomer_b_sp3_zmatrix(self):

      zmatrix = '''\
          H21    :1    DISTANCE     :2  ANGLE        :3 DIHEDRAL
          X21   H21   1.0000        :1  90.0000        :2    0.0000
          C21   H21  1.1066  X21  90.0000        :2  180.0000
          C22   C21  1.5353  H21  110.8679        :1  180.0000
          C23   C22  1.5178   C21  102.4718   H21 -149.3354
          C24   C23  1.3546   C22  111.3707   C21   16.4054
          C25   C21  1.5505   C22  105.0372   C23  -25.6092
          H22   C24  1.1007   C23  124.7945   C22  176.9201
          H23   C25  1.1007   C21  111.3707   C22  -92.1016
          H24   C25  1.1007   C21  112.5420   C22  147.5060
          H25   C23  1.1007   C22  123.7585   C21 -166.6366
          H26   C22  1.1007   C23  113.1858   C24  137.8623
          H27   C22  1.1007   C23  110.0838   C24 -102.2965
          H28   C21  1.1007   C22  108.7004   C23   90.5856
          X21   H21  1.0000   C21   90.0000   C22  180.0000
          0 1
        '''

      atom_name = [
        'H41', 'C4', 'C3', 'C2', 'C1', 'C5', 'H1', 'H51', 'H52', 'H2', 'H31', 'H32', 'H42'
      ]

      return textwrap.dedent(zmatrix), atom_name



