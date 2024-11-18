#!/usr/bin/env python3
#
# Lennard-Jones-Drill-2: Piperidine
# ----------------------------------

# Imports
# -------

import textwrap

class Piperidine(object):

  def __init__(self):

    self.resi_name = 'PIP'

  def get_monomer_a_species(self):

      '''

      Get the Monomer A Species

      '''

      monomer_a_species = {
          'H1': self.get_monomer_a_nitrogen_zmatrix(),
          'H2': self.get_monomer_a_carbon_zmatrix()
      }

      return monomer_a_species

  def get_monomer_a_nitrogen_zmatrix(self):

      zmatrix='''\
          H11
          N11   H11   0.9997
          C11   N11   1.4512   H11  111.1879
          C12   C11   1.5387   N11  108.9460    H11  164.2624
          C13   C12   1.5398   C11  110.8586    N11   30.0101
          C14   C13   1.5287   C12  110.2851    C11   31.8471
          C15   C14   1.5403   C13  110.6647    C12  -61.8944
          H12   C13   1.0858   C12  109.6180    C11  -88.2686
          H13   C13   1.0858   C12  110.3325    C11  154.1912
          H14   C14   1.0858   C13  111.0800    C12  176.0116
          H15   C14   1.0858   C13  108.9832    C12   59.1316
          H16   C15   1.0858   N11  107.9166    C11  160.1408
          H17   C15   1.0923   N11  112.2776    C11  -82.1438
          H18   C12   1.0858   C11  109.8805    N11  -91.1206
          H19   C12   1.0858   C11  109.3368    N11  151.8410
          H20   C11   1.0937   C12  110.0403    C13  -93.3694
          H21   C11   1.0858   C12  110.3187    C13  149.2429
          X11   H11   1.0000   N11   90.0000    C11    0.0000
          0 1
      '''

      atom_name = [
        'H11', 'N1', 'C2', 'C3', 'C4', 'C5', 'C6', 'H41', 'H42', 'H51', 'H52', 'H51', 'H52', 'H31', 'H32', 'H21', 'H22'
      ]

      return textwrap.dedent(zmatrix), atom_name

  def get_monomer_a_carbon_zmatrix(self):

      zmatrix = '''\
          H11
          C11   H11   1.0839
          N11   C11   1.4512   H11  111.1879
          C12   C11   1.5387   N11  108.9460    H11  164.2624
          C13   C12   1.5398   C11  110.8586    N11   30.0101
          C14   C13   1.5287   C12  110.2851    C11   31.8471
          C15   C14   1.5403   C13  110.6647    C12  -61.8944
          H12   C13   1.0848   C12  109.6180    C11  -88.2686
          H13   C13   1.0858   C12  110.3325    C11  154.1912
          H14   C14   1.0853   C13  111.0800    C12  176.0116
          H15   C14   1.0872   C13  108.9832    C12   59.1316
          H16   C15   1.0854   N11  107.9166    C11  160.1408
          H17   C15   1.0923   N11  112.2776    C11  -82.1438
          H18   C12   1.0851   C11  109.8805    N11  -91.1206
          H19   C12   1.0868   C11  109.3368    N11  151.8410
          H20   C11   1.0937   C12  110.0403    C13  -93.3694
          H21   C11   1.0839   C12  110.3187    C13  149.2429
          H21   N11   0.9997   C11  111.1879    C12  164.2624
          X11   H11   1.0000   C11   90.0000    N11    0.0000
          0 1
      '''

      atom_name = [
        'H22', 'C2', 'N1', 'C3', 'C4', 'C5', 'C6', 'H41', 'H42', 'H51', 'H52', 'H51', 'H52', 'H31', 'H32', 'H21', 'H11',
      ]

      return textwrap.dedent(zmatrix), atom_name

  def get_monomer_b_species(self):

      monomer_b_species = {
      }

      return monomer_b_species

