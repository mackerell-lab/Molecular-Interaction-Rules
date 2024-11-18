#!/usr/bin/env python3
#
# Lennard-Jones-Drill-2: Indole
# -----------------------------

# Imports
# -------
import textwrap

class Indole(object):

    def __init__(self):

        self.resi_name = 'INDO'

    def get_monomer_a_species(self):

        monomer_a_species = {
            'H1': self.monomer_a_nitrogen_hydrogen_zmatrix(),
            'RC1': self.monomer_a_aromatic_zmatrix()
        }

        return monomer_a_species

    def get_monomer_b_species(self):

        monomer_b_species = {
          'RC1': self.monomer_b_pi_stack_zmatrix(),
          'H1': self.monomer_b_nitrogen_hydrogen_zmatrix()
        }

        return monomer_b_species

    def monomer_a_aromatic_zmatrix(self):

        zmatrix = '''\
            X11
            N11  X11  1.1020
            C11  N11  1.3836  X11   59.0000
            C12  C11  1.4104  N11  130.4152  X11  180.0000
            C13  C12  1.4006  C11  117.2345  N11 -180.0000
            C14  C13  1.4233  C12  121.2757  C11   -0.0000
            C15  C14  1.3990  C13  121.3172  C12    0.0000
            C16  C11  1.4341  C12  122.4965  C13   -0.0000
            C17  C16  1.4395  C11  107.0716  C12 -180.0000
            C18  C17  1.3901  C16  106.8762  C15  180.0000
            H11  N11  1.0122  C11  125.2840  C12  0.0000
            H12  C14  1.0940  C13  119.0452  C12  180.0000
            H13  C15  1.0943  C14  120.6975  C13  180.0000
            H14  C17  1.0889  C16  127.5288  C11  180.0000
            H15  C18  1.0890  C17  130.1570  C16 -180.0000
            H16  C13  1.0939  C12  119.3067  C11  180.0000
            H17  C12  1.0944  C13  121.1026  C14 -180.0000
            0 1
        '''

        atom_name = [
          'NE1', 'CE2', 'CZ2', 'CH2', 'CZ3', 'CE3', 'CD2', 'CG', 'CD1',
          'HE1', 'HZ2', 'HH2', 'HZ3', 'HE3', 'HG', 'HD1'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def monomer_a_nitrogen_hydrogen_zmatrix(self):

        zmatrix = '''\
            H11
            N11  H11  1.0122
            C11  N11  1.3836  H11  125.2840
            C12  C11  1.4104  N11  130.4152  H11    0.0000
            C13  C12  1.4006  C11  117.2345  N11 -180.0000
            C14  C13  1.4233  C12  121.2757  C11   -0.0000
            C15  C14  1.3990  C13  121.3172  C12    0.0000
            C16  C11  1.4341  C12  122.4965  C13   -0.0000
            C17  C16  1.4395  C11  107.0716  C12 -180.0000
            C18  C17  1.3901  C16  106.8762  C15  180.0000
            H12  C14  1.0940  C13  119.0452  C12  180.0000
            H13  C15  1.0943  C14  120.6975  C13  180.0000
            H14  C17  1.0889  C16  127.5288  C11  180.0000
            H15  C18  1.0890  C17  130.1570  C16 -180.0000
            H16  C13  1.0939  C12  119.3067  C11  180.0000
            H17  C12  1.0944  C13  121.1026  C14 -180.0000
            X11  H11  1.0000  N11   90.0000  C11  180.0000
            0 1
        '''

        atom_name = [
          'HE1', 'NE1', 'CE2', 'CZ2', 'CH2', 'CZ3', 'CE3', 'CD2', 'CG', 'CD1',
          'HZ2', 'HH2', 'HZ3', 'HE3', 'HG', 'HD1'
        ]

        return textwrap.dedent(zmatrix), atom_name

    def monomer_b_pi_stack_zmatrix(self):

        zmatrix = '''\
            X21   :1  DISTANCE  :2    ANGLE    :3   90.0000
            N21  X21  1.1020    :1   90.0000   :2    0.0000
            C21  N21  1.3836  X21   59.0000   :1   DIHEDRAL
            C22  C21  1.4104  N21  130.4152  X21  180.0000
            C23  C22  1.4006  C21  117.2345  N21 -180.0000
            C24  C23  1.4233  C22  121.2757  C21   -0.0000
            C25  C24  1.3990  C23  121.3172  C22    0.0000
            C26  C21  1.4341  C22  122.4965  C23   -0.0000
            C27  C26  1.4395  C21  107.0716  C22 -180.0000
            C28  C27  1.3901  C26  106.8762  C25  180.0000
            H22  C24  1.0940  C23  119.0452  C22  180.0000
            H23  C25  1.0943  C24  120.6975  C23  180.0000
            H24  C27  1.0889  C26  127.5288  C21  180.0000
            H25  C28  1.0890  C27  130.1570  C26 -180.0000
            H26  C23  1.0939  C22  119.3067  C21  180.0000
            H27  C22  1.0944  C23  121.1026  C24 -180.0000
            0 1
        '''

        atom_name = [
          'NE1', 'CE2', 'CZ2', 'CH2', 'CZ3', 'CE3', 'CD2', 'CG', 'CD1',
          'HE1', 'HZ2', 'HH2', 'HZ3', 'HE3', 'HG', 'HD1'
        ],

        return textwrap.dedent(zmatrix), atom_name

    def monomer_b_nitrogen_hydrogen_zmatrix(self):

      zmatrix = '''\
            H21   :1  DISTANCE  :2   ANGLE    :3   90.0000
            X21  H21  1.0000   :1   90.0000   :2    0.0000
            N21  H21  1.0122  X21   90.0000   :1  DIHEDRAL
            C21  N21  1.3836  H21  125.2840   :1    0.0000
            C22  C21  1.4104  N21  130.4152  H21    0.0000
            C23  C22  1.4006  C21  117.2345  N21 -180.0000
            C24  C23  1.4233  C22  121.2757  C21   -0.0000
            C25  C24  1.3990  C23  121.3172  C22    0.0000
            C26  C21  1.4341  C22  122.4965  C23   -0.0000
            C27  C26  1.4395  C21  107.0716  C22 -180.0000
            C28  C27  1.3901  C26  106.8762  C25  180.0000
            H22  C24  1.0940  C23  119.0452  C22  180.0000
            H23  C25  1.0943  C24  120.6975  C23  180.0000
            H24  C27  1.0889  C26  127.5288  C21  180.0000
            H25  C28  1.0890  C27  130.1570  C26 -180.0000
            H26  C23  1.0939  C22  119.3067  C21  180.0000
            H27  C22  1.0944  C23  121.1026  C24 -180.0000
            0 1
        '''

      atom_name = [
        'HE1', 'NE1', 'CE2', 'CZ2', 'CH2', 'CZ3', 'CE3', 'CD2', 'CG', 'CD1',
        'HZ2', 'HH2', 'HZ3', 'HE3', 'HG', 'HD1'
      ]

      return textwrap.dedent(zmatrix), atom_name

