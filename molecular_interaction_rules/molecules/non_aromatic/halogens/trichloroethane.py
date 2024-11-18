#!/usr/bin/env python3
#
# Lennard-Jones-Drill-2: Trifluoroethane
# --------------------------------------

# Imports
# -------

import textwrap

class Trichloroethane(object):

  def __init__(self):

    self.resi_name = 'TFET'

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
        C12 CL11 1.7853
        C11 C12  1.5227 CL11 109.7299
        CL12 C12 1.7853 C11 109.7299 CL11 120.0000
        CL13 C12 1.7853 C11 109.7299 CL11 -120.0000
        H11 C11 1.0964 C12 111.4558 CL11 -60.0000
        H12 C11 1.0964 C12 111.4558 CL11 180.0000
        H13 C11 1.0964 C12 111.4557 CL11 60.0001
        X11 CL11 1.0000 C12  90.0000 C11 180.0000
        0 1
      '''

      atom_name = [
        'CL11', 'C1', 'C2', 'CL12', 'CL13', 'H21', 'H22', 'H23'
      ]

      return textwrap.dedent(zmatrix), atom_name

  def get_monomer_b_species(self):

      monomer_b_species = {
        'CL1': self.get_monomer_b_chloro_zmatrix()
      }

      return monomer_b_species

  def get_monomer_b_chloro_zmatrix(self):

      zmatrix = '''\
          CL21  :1 DISTANCE :2  ANGLE   :3  DIHEDRAL
          X21 CL21 1.0000   :1  90.0000   :2   0.0000
          C22 CL21 1.7853  X21  90.0000   :1  180.0000
          C21 C22  1.5227 CL21 109.7299   :1  180.0000
          CL22 C22 1.7853 C21 109.7299 CL21 120.0000
          CL23 C22 1.7853 C21 109.7299 CL21 -120.0000
          H21 C21 1.0964 C22 111.4558 CL21 -60.0000
          H22 C21 1.0964 C22 111.4558 CL21 180.0000
          H23 C21 1.0964 C22 111.4557 CL21 60.0001
          0 1
        '''

      atom_name = [
        'CL11', 'C1', 'C2', 'CL12', 'CL13', 'H21', 'H22', 'H23'
      ]

      return textwrap.dedent(zmatrix), atom_name

