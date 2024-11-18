#!/usr/bin/env python3
#
# Lennard-Jones-Drill-2: Chloroethane
# -----------------------------------

# Imports
# -------

import textwrap

class Chloroethane(object):

  def __init__(self):

    self.resi_name = 'CLET'

  def get_monomer_a_species(self):

    '''

    Get the Monomer A Species

    '''

    monomer_a_species = {
      'CL1': self.get_monomer_a_chloro_zmatrix()
    }

    return monomer_a_species

  def get_monomer_a_chloro_zmatrix(self):

      zmatrix = '''\
          CL11
          C12 CL11 1.7767
          C11 C12  1.5143 CL11 110.5148
          H11 C12  1.0936 C11  110.7942  CL11 119.6436
          H12 C12  1.0936 C11 110.7942 H11 120.7128
          H13 C11  1.0946 C12 110.1540 H11 60.3561
          H14 C11  1.0950 C12 111.2308 H11 180.0000
          H15 C11  1.0950 C12 111.2308 H11 -59.2502
          X11 CL11 1.0000 C12 90.0000  C11 180.0000
          0 1
      '''

      atom_name = [
        'CL11', 'C1', 'C2', 'H11', 'H12', 'H21', 'H22', 'H23'
      ]

      return textwrap.dedent(zmatrix), atom_name

  def get_monomer_b_species(self):

      monomer_b_species = {
        'CL1': self.get_monomer_b_chloro_zmatrix()
      }

      return monomer_b_species

  def get_monomer_b_chloro_zmatrix(self):

      zmatrix = '''\
            CL21  :1 DISTANCE  :2 ANGLE  :3  DIHEDRAL
            X21 CL21 1.0000  :1  90.0000    :2   0.0000
            C22 CL21 1.7767  X21  90.0000   :1  180.0000
            C21 C22  1.5143 CL21 110.5148   :1  180.0000
            H21 C22  1.0936 C21  110.7942  CL21 119.6436
            H22 C22  1.0936 C21 110.7942 H21 120.7128
            H23 C21  1.0946 C22 110.1540 H21 60.3561
            H24 C21  1.0950 C22 111.2308 H21 180.0000
            H25 C21  1.0950 C22 111.2308 H21 -59.2502
            0 1
        '''

      atom_name = [
        'CL11', 'C1', 'C2', 'H11', 'H12', 'H21', 'H22', 'H23'
      ]

      return textwrap.dedent(zmatrix), atom_name
