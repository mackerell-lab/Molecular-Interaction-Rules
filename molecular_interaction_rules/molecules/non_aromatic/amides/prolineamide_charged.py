#!/usr/bin/env python3
#
# Lennard-Jones-Drill-2: ProlineamideCharged
# ------------------------------------------

# Imports
# -------

import textwrap

class ProlineamideCharged(object):

  '''

  RESI PNH2          1.00 ! C5H11N2O, prolineamide (aka TP2, PAMD), R. Dunbrack
  ! charges adjusted for consistency with CGenFF charges
  ATOM N    NG3P2   -0.22
  ATOM HN1  HGP2     0.34
  ATOM HN2  HGP2     0.34
  ATOM CD   CG3C54  -0.35
  ATOM CB   CG3C52  -0.12 !    HN1   HD1 HD2
  ATOM CG   CG3C52  -0.12 !      \    \ /
  ATOM CA   CG3C53   0.12 !  HN2--N---CD   HG1
  ATOM C    CG2O1    0.51 !       |(+)  \  /
  ATOM O    OG2D1   -0.51 !       |      CG
  ATOM HA   HGA1     0.09 !       |     / \
  ATOM HB1  HGA2     0.09 !    HA-CA--CB   HG2
  ATOM HB2  HGA2     0.09 !       |   / \
  ATOM HG1  HGA2     0.09 !       | HB1 HB2
  ATOM HG2  HGA2     0.09 !     O=C
  ATOM HD1  HGA2     0.28 !       |
  ATOM HD2  HGA2     0.28 !       NT
  ATOM NT   NG2S2   -0.62 !      / \
  ATOM HT1  HGP1     0.31 !   HT1   HT2
  ATOM HT2  HGP1     0.31

  Rule 3

  '''

  __CGENFF_ATOM_TYPES__ = {
    'H1': ['HGA1', 'CG3C53'],
    # 'H2': ['HGA2', 'CG3C54']
  }

  __DGENFF_ATOM_TYPES__ = {
  }

  def __init__(self):

    self.resi_name = 'PNH2'

  def get_monomer_a_species(self):

    '''

    Get the Monomer A Species

    '''

    monomer_a_species = {
        'H1': self.get_monomer_a_hydrogen_zmatrix()
    }

    return monomer_a_species

  def get_monomer_a_hydrogen_zmatrix(self):

    zmatrix = '''\
      H11
      C12  H11 1.0846
      C11  C12 1.5300   H11  104.5232
      N11  C11 1.3445   C12  115.3149  H11  -97.9238
      N12  C12 1.4641   C11  112.2071  N11   20.4696
      H12  N12 1.0010   C12  109.4950  C11 -143.0365
      H13  N12 1.0010   C12  109.4950  C11   37.0365
      C13  N12 1.4615   C12  105.8259  C11   99.2627
      H14  C13 1.0861   N12  108.5988  C12  -81.7117
      H15  C13 1.0830   N12  110.5262  C12  160.4969
      C14  C13 1.5322   N12  105.5222  H12  -81.1190
      H16  C14 1.0841   C13  112.8941  N12 -156.2638
      H17  C14 1.0864   C13  109.7458  N12   83.8645
      C15  C12 1.5538   C11  111.9567  N11  141.5078
      H18  C15 1.0807   C12  110.2896  C11   -0.3212
      H19  C15 1.0828   C12  111.1976  C11  118.1844
      O11  C11 1.2019   C12  121.1157  N12 -162.3860
      H20  N11 0.9942   C11  118.5945  C12 -179.6150
      H21  N11 0.9935   C11  120.0679  C12   -8.4508
      X11  H11 1.0000   C12   90.0000  C11    0.0000
      1 1
    '''

    atom_name = [
      'HA', 'CA', 'C', 'NT', 'N', 'HN1', 'HN2', 'CD', 'HD1', 'HD2', 'CG', 'HG1', 'HG2', 'CB', 'HB1', 'HB2', 'O', 'HT1', 'HT2'
    ]

    return textwrap.dedent(zmatrix), atom_name

  def get_monomer_b_species(self):

    '''

    Get the Monomer B Species

    '''

    monomer_b_species = {
      'H1': self.get_monomer_b_hydrogen_zmatrix()
    }

    return monomer_b_species

  def get_monomer_b_hydrogen_zmatrix(self):

    zmatrix = '''\
      H21   :1  DISTANCE  :2 90.0000    :3   90.0000
      X21  H21 1.0000     :1 90.0000    :2    0.0000
      C22  H21 1.0846   X21  90.0000    :2  180.0000
      C21  C22 1.5300   H21  104.5232   :1  180.0000
      N21  C21 1.3445   C22  115.3149  H21  -97.9238
      N22  C22 1.4641   C21  112.2071  N21   20.4696
      H22  N22 1.0010   C22  109.4950  C21 -143.0365
      H23  N22 1.0010   C22  109.4950  C21   37.0365
      C23  N22 1.4615   C22  105.8259  C21   99.2627
      H24  C23 1.0861   N22  108.5988  C22  -81.7117
      H25  C23 1.0830   N22  110.5262  C22  160.4969
      C24  C23 1.5322   N22  105.5222  H22  -81.1190
      H26  C24 1.0841   C23  112.8941  N22 -156.2638
      H27  C24 1.0864   C23  109.7458  N22   83.8645
      C25  C22 1.5538   C21  111.9567  N21  141.5078
      H28  C25 1.0807   C22  110.2896  C21   -0.3212
      H29  C25 1.0828   C22  111.1976  C21  118.1844
      O21  C21 1.2019   C22  121.1157  N22 -162.3860
      H30  N21 0.9942   C21  118.5945  C22 -179.6150
      H31  N21 0.9935   C21  120.0679  C22   -8.4508
      1 1
    '''

    atom_name = [
      'HA', 'CA', 'C', 'NT', 'N', 'HN1', 'HN2', 'CD', 'HD1', 'HD2', 'CG', 'HG1', 'HG2', 'CB', 'HB1', 'HB2', 'O', 'HT1', 'HT2'
    ]

    return textwrap.dedent(zmatrix), atom_name

