#!/usr/bin/env python3
#
# Lennard-Jones-Drill-2: Pyridine
# -------------------------------

# Imports
# -------
import textwrap

class Pyridinium(object):

  def __init__(self):

    self.resi_name = 'PIUM'

  def get_monomer_a_species(self):

    '''

    Get the Monomer A Species

    '''

    monomer_a_species = {
      'H1': self.monomer_a_hydrogen_zmatrix()
    }

    return monomer_a_species

  def get_monomer_b_species(self):

    monomer_b_species = {
      'H1': self.monomer_b_hydrogen_zmatrix()
    }

    return monomer_b_species

  def monomer_a_hydrogen_zmatrix(self):

      zmatrix = '''\
           H11
           N11  H11  0.9996
           C11  N11  1.3208  H11  180.0000
           C12  C11  1.3774  N11  124.0000  H11   0.0000
           C13  C12  1.3774  C11  120.0000  N11   0.0000
           C14  C13  1.3774  C12  120.0000  C11   0.0000
           C15  N11  1.3208  C11  116.0000  C12   0.0000
           H11  C11  1.0756  C12  120.0000  C13 180.0000
           H12  C12  1.0756  C13  120.0000  C14 180.0000
           H13  C13  1.0756  C14  120.0000  C15 180.0000
           H14  C14  1.0756  C15  120.0000  N11 180.0000
           H15  C15  1.0756  N11  120.0000  C11 180.0000
           X11  H11  1.0000  N11   90.0000  C11 180.0000
           1 1
        '''

      atom_name = [
         'H6', 'N1', 'C1', 'C2', 'C3', 'C4', 'C5', 'H1', 'H2', 'H3', 'H4', 'H5'
      ]

      return textwrap.dedent(zmatrix), atom_name

  def monomer_b_hydrogen_zmatrix(self):

    zmatrix = '''\
           H21   :1  DISTANCE   :2   ANGLE   :3   90.0000
           X21   H21  1.0000    :1   90.0000    :2   DIHEDRAL
           N21  H21  0.9996  X21   90.0000      :1  180.0000
           C21  N21  1.3208  H21  180.0000   X21     0.0000
           C22  C21  1.3774  N21  124.0000  H21   0.0000
           C23  C22  1.3774  C21  120.0000  N21   0.0000
           C24  C23  1.3774  C22  120.0000  C21   0.0000
           C25  N21  1.3208  C21  116.0000  C22   0.0000
           H21  C21  1.0756  C22  120.0000  C23 180.0000
           H22  C22  1.0756  C23  120.0000  C24 180.0000
           H23  C23  1.0756  C24  120.0000  C25 180.0000
           H24  C24  1.0756  C25  120.0000  N21 180.0000
           H25  C25  1.0756  N21  120.0000  C21 180.0000
           1 1
    '''

    atom_name = [
      'H6', 'N1', 'C1', 'C2', 'C3', 'C4', 'C5', 'H1', 'H2', 'H3', 'H4', 'H5'
    ]

    return textwrap.dedent(zmatrix), atom_name


